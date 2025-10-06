#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import dendropy
from scipy.optimize import minimize

from utils import (
    tree_llh_continuous,        # llh object factory
    calc_llh_per_tree,           # expects (rates, llhobj) -> float
    read_nwk_trees,             # returns list[str] newick
    nwk2tree,                   # newick str -> dendropy.Tree
    parse_trees,                # path -> dict[newick_str -> weight]
    segment_mapping             # annotate gene tree edges with sp_segment
)

# ---------- Objective helpers ----------
def objective_sum_llh(rates: np.ndarray, llhobjs: List, pool: Pool) -> float:
    vals = pool.map(calc_llh_per_tree, [(rates, obj) for obj in llhobjs])
    return float(np.sum(vals))


def build_llhobjs_and_init(
    tip_trait: Dict[str, float],
    trees: List,
) -> Tuple[List, float, float]:
    """
    Build llh objects, initialize pruning, and compute initial values.
    """
    llhobjs = [tree_llh_continuous(t) for t in trees]

    # Max root distance per tree; then global max height
    tree_heights = []
    for llh in llhobjs:
        dists = llh.tree.calc_node_root_distances()
        tree_heights.append(float(np.max(list(dists.values()) if isinstance(dists, dict) else dists)))
    max_height = float(np.max(tree_heights))

    # Initialize pruning with the global max height
    for obj in llhobjs:
        obj.initialize_pruning(tip_trait, max_height)

    # Initial guess
    init_value = float(np.var(list(tip_trait.values())) / (max_height if max_height > 0 else 1.0))
    return llhobjs, max_height, init_value


def infer_rates_for_trait(
    tip_trait: Dict[str, float],
    trees: List[dendropy.Tree],
    segment2rate: Dict[int, int],
    trait_idx: int,
    condition: str,
    num_process: int
) -> None:

    # Build llh objects + initialization
    llhobjs, _, init_scalar = build_llhobjs_and_init(tip_trait, trees)

    # Number of distinct rate families (1-indexed in segment2rate; we shift to 0-index elsewhere)
    numrate = int(len(np.unique(list(segment2rate.values()))))

    # Initial vector
    x0 = np.full(shape=(numrate,), fill_value=max(init_scalar, 1e-6), dtype=np.float64)

    # Optimizer + parallel evaluation
    pool = Pool(num_process)
    try:
        pool = Pool(num_process)
        obj = partial(objective_sum_llh, llhobjs=llhobjs, pool=pool)
        res = minimize(
            obj,
            x0,
            method="Nelder-Mead",
            bounds=((1e-8, np.inf),) * numrate,
        )
    finally:
        pool.close()
        pool.join()

    record = {
    "n_trees": len(trees),
    "minimize_success": bool(res.success),
    **{f"rate_{i}": float(v) for i, v in enumerate(res.x)},
    "condition": condition,
    "trait_idx": int(trait_idx),
    }

    return record


def main():
    parser = argparse.ArgumentParser(description="Rate inference over gene trees.")
    parser.add_argument("-f", "--traits", help="CSV of traits (rows = traits, cols = taxa)", type=Path, required=True)
    parser.add_argument("-g", "--genetrees", help="gene trees (Newick) path; branch lengths in 2N units", type=Path, required=True)
    parser.add_argument("-s", "--speciestree", help="species tree (Newick) path; branch lengths in 2N units", type=Path, required=True)
    parser.add_argument("-r", "--ratetree", help="rate-family tree (same topology as species tree)", type=Path, required=True)
    parser.add_argument("-c", "--condition", help="label for outputs", type=str, default="")
    parser.add_argument("-n", "--num_process", help="processes for parallel objective eval", type=int, default=1)
    

    args = parser.parse_args()

    # Basic validation
    for p in [args.traits, args.genetrees, args.speciestree, args.ratetree]:
        if p is None or not Path(p).exists():
            raise FileNotFoundError(f"Missing or unreadable path: {p}")
        
    # Species tree 
    sptree = nwk2tree(read_nwk_trees(str(args.speciestree))[0])

    speciation_time = sorted(
        {node.distance_from_tip() for node in sptree.nodes() if len(node.child_nodes()) == 2},
        reverse=True
    )

    # Rate-family mapping (species-tree edges -> family id [1..K])
    ratetree = dendropy.Tree.get(
        data=read_nwk_trees(str(args.ratetree))[0],
        schema="newick",
        taxon_namespace=sptree.taxon_namespace
    )

    for idx, edge in enumerate(sptree.preorder_edge_iter()):
        setattr(edge, "sp_segment", idx)

    segment2rate = {}
    for sedge, redge in zip(sptree.preorder_edge_iter(), ratetree.preorder_edge_iter()):
        # Store family id as integer (expects branch length encodes family id)
        segment2rate[getattr(sedge, "sp_segment")] = int(round(redge.length))

    # Gene trees with weights
    tree_weights = parse_trees(str(args.genetrees))  # dict[newick -> weight]
    genetrees: List[dendropy.Tree] = []
    for nwk, freq in tree_weights.items():
        t = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=sptree.taxon_namespace)
        setattr(t, "weight", float(freq))
        genetrees.append(t)

    # Assign rate families to segments in each gene tree
    for t in genetrees:
        segment_mapping(sptree, t, speciation_time)
        for node in t.nodes():
            if hasattr(node.edge, "sp_segment"):
                # convert 1..K -> 0..K-1 for indexing elsewhere if needed
                setattr(node, "rate_family", int(segment2rate[node.edge.sp_segment]) - 1)

    # Traits
    traits_df = pd.read_csv(args.traits, dtype=float)
    tip_traits = [traits_df.loc[i].to_dict() for i in range(traits_df.shape[0])]
    
    # Output
    out_path = f"./{args.condition}.json"
    records = []

    # Run inference per trait
    for idx, tip_trait in enumerate(tip_traits):
        record = infer_rates_for_trait(
                    tip_trait=tip_trait,
                    trees=genetrees,
                    segment2rate=segment2rate,
                    trait_idx=idx,
                    condition=args.condition,
                    num_process=args.num_process,
                )
        records.append(record)
        print(record)

    with open(out_path, "w") as f:
        for item in records:
            json_line = json.dumps(item)
            f.write(json_line + "\n")

if __name__ == "__main__":
    main()
