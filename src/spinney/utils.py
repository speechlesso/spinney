import numpy as np
import pandas as pd
from itertools import combinations
import copy, json, re
import dendropy

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import seaborn as sns


# %% TREE process
def read_nwk_trees(treefile):
    # read newick trees from a file
    nwks = []
    with open(treefile, "r") as f:
        for line in f.readlines():
            if line != "\n" and len(re.findall("\(", line)) != 0:
                nwks.append(line.split("\n")[0])
    return nwks


def nwk2tree(newick_tree, is_4N_unit=False):
    tree = dendropy.Tree.get(data=newick_tree, schema="newick")
    # numInternalNode = 1
    internal_node = [x for x in tree.postorder_internal_node_iter()]
    for clade in internal_node:
        if clade.edge_length is None:
            clade.edge_length = 0.0
    if is_4N_unit:
        for x in tree.postorder_node_iter():
            x.edge_length *= 2
    return tree


# Add rate as attribute for each branch segment with a rate tree
def branch_rate_family(treeobj, sptree, rate_tree=None):
    if rate_tree != None:  # multiple rate
        # rate_tree has identical topology as speceis tere
        [
            setattr(sedge, "rate_family", int(redge.length))
            for sedge, redge in zip(
                [se for se in sptree.preorder_edge_iter()],
                [re for re in rate_tree.preorder_edge_iter()],
            )
        ]
    else:
        # single rate
        [setattr(clade, "rate_family", 0) for clade in treeobj.nodes()]
    return treeobj


def parse_trees(treefile):
    with open(treefile, "r") as f:
        lines = f.readlines()
    # Remove any whitespace or newline characters from each line
    lines = [line.strip() for line in lines]

    # Check where the gene trees and genetree freqs sections start and end
    begin_genetrees = lines.index("begin genetrees;")
    end_genetrees = lines.index("end;")
    begin_genetreefreqs = lines.index("begin genetreefreqs;")

    # Extract gene trees and their frequencies
    genetreesnwk = lines[begin_genetrees + 1 : end_genetrees]
    genetreefreqs = lines[begin_genetreefreqs + 1 : -1]  # -1 to exclude the "end;" line

    # Convert frequencies from string to float
    genetreefreqs = [float(freq) for freq in genetreefreqs]

    tree_weights = dict(zip(genetreesnwk, genetreefreqs))

    return tree_weights


# Add speciation time slice as internal nodes to gene trees
def time_slice_node(tree, slice2tip):
    root2tip = tree.max_distance_from_root()
    root2slice = root2tip - slice2tip
    if root2slice < 0 or slice2tip < 0:
        print(
            "distance for root to time slice is %.2f, distance from time slice to tip is %.2f"
            % (root2slice, slice2tip)
        )
        return tree

    # get possible parents for time-slice node in each lineage
    parent4slice = [
        x for x in tree.levelorder_node_iter() if x.distance_from_tip() >= slice2tip
    ]

    # add time slice as new nodes
    for parent_node in parent4slice:
        slice2parent = parent_node.distance_from_tip() - slice2tip
        children = parent_node.child_nodes()
        if len(children) > 1:
            for child in children:
                new_edge = child.edge_length - slice2parent
                if new_edge >= 0:
                    child = parent_node.remove_child(child)
                    child.edge_length = new_edge
                    slice_node = parent_node.new_child(edge_length=slice2parent)
                    slice_node.add_child(child)
                    setattr(slice_node, "time_slice", True)
                    tree.update_bipartitions(suppress_unifurcations=False)

    for idx, x in enumerate(tree.preorder_internal_node_iter()):
        x.label = "N%s" % idx
        try:
            x.time_slice
        except AttributeError:
            setattr(x, "time_slice", False)

    for idx, x in enumerate(tree.leaf_node_iter()):
        x.label = x.taxon.label
        setattr(x, "time_slice", True)

    return tree


def find_path_to_root(target):
    path = [target.edge]
    node = target
    while node.parent_node is not None:
        node = node.parent_node
        path.append(node.edge)
    return path


def segment_mapping(sptree, genetree, speciation_time):
    """
    Assign rate from the species tree to each segment of a gene tree
    """
    # speciation_time.sort()
    # input species tree, map each edge to an integer
    leaf_names = [x.taxon.label for x in sptree.leaf_nodes()]
    leaf_names.sort()
    sedges = [edge for edge in sptree.preorder_edge_iter()]
    sedge2int = {edge: idx for idx, edge in enumerate(sedges)}

    # get path from root to edge, each element is an edge object in dendropy
    gpaths = [find_path_to_root(leaf) for leaf in genetree.leaf_node_iter()]

    for gpath in gpaths:
        leaf_taxon = gpath[0].head_node.taxon.label
        spath = [
            find_path_to_root(leaf)
            for leaf in sptree.leaf_node_iter()
            if leaf.taxon.label == leaf_taxon
        ][0]
        # edge x : root------ tail ---x--- head
        for edge in gpath:
            try:
                # avoid redundant assigning
                edge.sp_segment
            except AttributeError:
                # find mapped branch in sptree
                idx = [edge.head_node.distance_from_tip() >= t for t in speciation_time]

                if any(idx):  # non-terminal edge
                    distance2tip = max(
                        [t for t, flag in zip(speciation_time, idx) if flag]
                    )
                else:  # terminal edge
                    distance2tip = 0.0

                mapped_sedge = [
                    se
                    for se in spath
                    if se.head_node.distance_from_tip() <= distance2tip
                ]
                setattr(edge, "sp_segment", sedge2int[mapped_sedge[-1]])

    return genetree


def common_ancestor(tree, targets):
    """Most recent common ancestor (clade) of all the given targets.
       Edge cases:
        - If no target is given, returns self.root
        - If 1 target is given, returns the target
        - If any target is not found in this tree, raises a ValueError
    Modified from Bio.Phylo.BaseTree.common_ancestor for dendropy tree object
    """
    paths = [[x.head_node for x in find_path_to_root(t)[:-1]] for t in targets]
    for p, t in zip(paths, targets):
        if p is None:
            raise ValueError(f"target {t!r} is not in this tree")
    mrca = tree.nodes()[0]

    [p.reverse() for p in paths]

    for level in zip(*paths):
        ref = level[0]
        for other in level[1:]:
            if ref is not other:
                break
        else:
            mrca = ref
        if ref is not mrca:
            break
    return mrca


# %% Class for calculating likelihood for continuous trait
class tree_llh_continuous:
    def __init__(self, treeobj, vec_size=100):
        self.tree = copy.deepcopy(treeobj)
        self.tree.root = self.tree.nodes()[0]
        self.leaves = sorted(self.tree.leaf_nodes(), key=lambda x: x.taxon.label)
        self.ordered_vertices = [x for x in self.tree.postorder_node_iter()]

        self.node_to_num = {
            n: i for i, n in enumerate([x for x in self.tree.postorder_node_iter()])
        }
        self.num_to_node = {v: k for k, v in self.node_to_num.items()}
        self.stationary = np.ones(vec_size) / vec_size

    def C_mat(self):
        # construct covariance matrix
        # diag elements
        leaves = sorted(self.tree.leaf_nodes(), key=lambda x: x.taxon.label)
        cov = np.diag(self.tree.calc_node_root_distances())
        # off-diag elements
        tip_to_num = {n.taxon: i for i, n in enumerate(leaves)}
        leaves_pair = [(t1, t2) for t1, t2 in combinations(leaves, 2)]

        for pair in leaves_pair:
            mrca = common_ancestor(self.tree, pair)
            r, c = [tip_to_num[x.taxon] for x in pair]
            cov[r, c] = cov[c, r] = mrca.distance_from_root()
        return cov

    def initialize_pruning(self, tip_trait, tree_height, vec_size=100):
        """
        INITIALIZE FOR ONE TRAIT EVERY TIME!
        THE TRAIT RANGE VARIES DEPENDING ON TIP TRAITS!

        tip_trait:a dict {leavename:single_trait_value}
        vec_size: vector length for range of trait
        """

        # discretize continous traits
        max_val = max(([v for v in tip_trait.values()]))
        min_val = min(([v for v in tip_trait.values()]))
        dev = np.sqrt(tree_height)
        self.trait_vec = np.linspace(
            min_val - 3 * dev, max_val + 3 * dev, vec_size
        ).flatten()

        # calculate state of tips
        tip_state = {}
        for key, val in tip_trait.items():
            temp_vec = np.zeros(self.trait_vec.size)
            # trait value matches one in discretized vector exactly
            val_index = np.argwhere(self.trait_vec == val).flatten()
            if len(val_index) != 0:
                temp_vec[val_index] = 1.0
            else:
                # trait value falls btween two numbers corresponding to trait range, e.g. 0.5 in [0,1]
                upperind = np.min(np.argwhere(self.trait_vec > val))
                lowerind = upperind - 1
                stepsize = self.trait_vec[upperind] - self.trait_vec[lowerind]
                weight = (val - self.trait_vec[lowerind]) / stepsize
                temp_vec[upperind] = weight
                temp_vec[upperind - 1] = 1 - temp_vec[upperind]
            tip_state.update({key: temp_vec})

        for leaf in self.tree.leaf_nodes():
            setattr(leaf, "state", tip_state[leaf.taxon.label])

    def calc_llh_with_tree(self, rate):
        nstate = self.stationary.size
        # initialize state for internal nodes
        for clade in self.tree.internal_nodes():
            setattr(clade, "state", np.zeros(nstate))

        # cache probability density matrix for trees
        for clade in self.ordered_vertices:
            if clade == self.tree.root:
                setattr(clade, "probmat", np.eye(nstate))
            else:
                prob_mat = np.zeros((nstate, nstate))
                for row in range(nstate):
                    for col in range(row + 1):
                        prob_mat[row, col] = bm_prob_density(
                            self.trait_vec[col],
                            self.trait_vec[row],
                            rate[clade.rate_family],
                            clade.edge_length,
                        )
                        prob_mat[col, row] = prob_mat[row, col]

                dx = self.trait_vec[1] - self.trait_vec[0]
                prob_mat *= dx
                setattr(clade, "probmat", prob_mat)

        for parent in self.tree.postorder_internal_node_iter():
            child_llh = []
            for i, child in enumerate(parent.child_nodes()):
                mat = child.probmat
                vec = child.state
                prod = mat @ vec
                child_llh.append(prod)

            if len(child_llh) == 1:
                parent.state = child_llh[0]
            else:
                product = child_llh[0] * child_llh[1]
                parent.state = product

        # warning raising if the rate is close to zero
        negative_llh = -np.log(self.stationary.dot(self.tree.root.state))
        return negative_llh


# %% HELPER FUNCTION FOR LIKLIHOOD
def bm_prob_density(x0, x1, rate, t):
    # x0: initial value; x1: final value; rate:sigma2; t:branch length
    var = max(rate * t, 1e-5)  # avoid extremely large prob density
    return (1.0 / np.sqrt(2 * np.pi * var)) * np.exp(-((x1 - x0) ** 2) / (2 * var))


def mle_sigma2(leaves, tip_trait, cov):
    leaves_name = [x.taxon.label for x in leaves]
    n_taxa = len(tip_trait)
    ones = np.ones([n_taxa, 1])
    x = np.array([tip_trait[x] for x in leaves_name]).reshape(n_taxa, 1)
    c = cov
    inv_c = np.linalg.inv(c)
    a = (np.linalg.inv(ones.T @ inv_c @ ones)) @ (ones.T @ inv_c @ x)
    return ((x - a).T @ inv_c @ (x - a) / n_taxa)[0][0]


def llh_surface(tip_trait, cov, sigma2):
    n_taxa = len(tip_trait)
    sort_key = [k for k in tip_trait.keys()]
    sort_key.sort()
    x = np.array([tip_trait[k] for k in sort_key]).reshape(n_taxa, 1)

    inv_c = np.linalg.inv(cov * sigma2)
    det_c = np.linalg.det(cov * sigma2)

    denom = np.exp(-0.5 * (x - np.mean(x)).T @ inv_c @ (x - np.mean(x)))
    nome = np.sqrt((2 * np.pi) ** n_taxa * det_c)
    return -np.log(denom / nome).flatten()[0]


def calc_llh_per_tree(x):
    rates, llhobj = x
    return llhobj.calc_llh_with_tree(rates) * llhobj.tree.weight


def init_llh_per_tree(x):
    llhobj, tip_trait, tree_height = x
    llhobj.initialize_pruning(tip_trait, tree_height)
    return llhobj


# %% continous trait simulator
def get_continous_trait(mean, cov_mat, seed=None):
    # construct cov_var matrix then sample trait for the tips.
    if seed != None:
        np.random.seed(seed)
    return np.random.multivariate_normal(mean, cov_mat)


# simulate continous trait with brownian motion
def sim_trait_val(species_tree, rate, root_trait, seed=None):
    leaves = species_tree.get_terminals()
    n_taxa = len(leaves)
    # diag elements
    diag = species_tree.distance(leaves[0])
    cov = np.diag([diag,]* n_taxa)

    # off-diag elements
    internal = [x for x in species_tree.get_nonterminals() if x != species_tree.root]
    tip_to_num = {n.name: i for i, n in enumerate(leaves)}

    leaves_pair = [x for x in combinations(leaves, 2)]
    for pair in leaves_pair:
        ca = species_tree.common_ancestor(pair)
        r, c = [tip_to_num[x.name] for x in pair]
        cov[r, c] = cov[c, r] = rate * species_tree.distance(ca)
    tip_init = dict(
        zip([x.name for x in leaves], get_continous_trait(root_trait, cov, seed))
    )  # {leafname:trai_value}

    return tip_init


def sim_ancestral_state(tree, rate, root_dist=None, seed=None):
    np.random.seed(seed)
    if root_dist == None:
        root_dist = 0

    for parent in tree.find_clades():
        if parent == tree.root:
            setattr(parent, "simtrait", root_dist)
        for child in parent.clades:
            trait_temp = np.random.normal(
                loc=parent.simtrait,
                scale=np.sqrt(child.branch_length * rate[child.rate_family]),
            )
            setattr(child, "simtrait", trait_temp)


# %% Ancestral states

# Path / I/O utilities
@dataclass
class Paths:
    folder: Path
    species_tree_file: Path
    rate_tree_file: Path
    trait_file: Path
    genetrees_file: Path
    rate_file: Path


def build_paths(prefix: str, treetype: str, base_dir: str ) -> Paths:
    base = Path(base_dir)
    species_tree_file = base / f"{prefix}.st.txt"
    rate_tree_file = base / f"{prefix}.ratetree_1p.nwk"
    trait_file = base / f"{prefix}.csv"
    genetrees_file = base / f"{prefix}.{treetype}.txt"
    rate_file = base / f"{prefix}_{treetype}_1p.json"
    return Paths(
        base, species_tree_file, rate_tree_file, trait_file, genetrees_file, rate_file
    )


def load_species_and_rate_trees(
    paths: Paths,
) -> Tuple[dendropy.Tree, dendropy.Tree, Dict[int, int]]:
    sptree = nwk2tree(read_nwk_trees(str(paths.species_tree_file))[0])
    ratetree = dendropy.Tree.get(
        data=read_nwk_trees(str(paths.rate_tree_file))[0],
        schema="newick",
        taxon_namespace=sptree.taxon_namespace,
    )
    # tag edges in species tree with integer segments
    for idx, edge in enumerate(sptree.preorder_edge_iter()):
        setattr(edge, "sp_segment", idx)
    segment2rate = {
        sedge.sp_segment: int(redge.length)
        for sedge, redge in zip(
            sptree.preorder_edge_iter(), ratetree.preorder_edge_iter()
        )
    }
    return sptree, ratetree, segment2rate


def load_gene_trees(paths: Paths, taxon_ns) -> List[dendropy.Tree]:
    # read with weights; parse_trees returns {nwk: weight}
    tree_weights = parse_trees(str(paths.genetrees_file))
    genetrees: List[dendropy.Tree] = []
    for nwk, freq in tree_weights.items():
        t = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=taxon_ns)
        t.weight = freq
        genetrees.append(t)
    return genetrees


def load_traits_and_rates(paths: Paths) -> Tuple[List[Dict], pd.DataFrame]:
    traits_df = pd.read_csv(paths.trait_file, comment="#", dtype=float)
    tip_traits = [traits_df.loc[i].to_dict() for i in range(traits_df.shape[0])]

    with open(paths.rate_file) as f:
        est = [json.loads(line) for line in f.readlines()]
    
    est = pd.DataFrame(est)
    est.drop_duplicates("trait_idx", inplace=True)
    est.sort_values(by="trait_idx", inplace=True)
    return tip_traits, est


# Gene-tree preparation
def slice_and_map_gene_trees(
    sptree: dendropy.Tree,
    genetrees: List[dendropy.Tree],
    segment2rate: Dict[int, int],
    time_grid: np.ndarray,
) -> None:
    for t in genetrees:
        curr_height = t.max_distance_from_root()
        for ts in time_grid[1:]:
            if ts < curr_height:
                time_slice_node(t, ts)
        segment_mapping(sptree, t, time_grid)
        for node in t.nodes():
            if hasattr(node.edge, "sp_segment"):
                setattr(node, "rate_family", segment2rate[node.edge.sp_segment] - 1)


def ancestral_state_map(tree, time_to_leaves):
    '''
    Given a tree and time slices, return a 2d dataframe
    whose row is time from tip, column is leaf names
    (i,j) entry represent ancestral node for leaf j at time i
    '''
    # Get all leaf nodes
    leaf_nodes = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter()}
    ancestor_table = {t: {leaf_name: None for leaf_name in leaf_nodes} for t in time_to_leaves}

    # Iterate through time slices and find ancestors

    for t in time_to_leaves:
        for leaf_name, leaf in leaf_nodes.items():
            node = leaf
            last_valid_ancestor = None  # Track the last valid ancestor before or at time t
            
            while node:
                curr_dist = node.distance_from_tip()
                # Move up while we haven't reached the time slice
                if (curr_dist <= t) or (curr_dist - t < 1e-6):  
                    last_valid_ancestor = node 
                    node = node.parent_node
                    
                else:
                    break
            
            ancestor_table[t][leaf_name] = last_valid_ancestor if last_valid_ancestor else None
            
            
            if last_valid_ancestor != tree.root:
                ancestor_table[t][leaf_name] = last_valid_ancestor if last_valid_ancestor else None
            # edge cases
            elif last_valid_ancestor.distance_from_tip() == t:
                ancestor_table[t][leaf_name] = last_valid_ancestor
            else:
                ancestor_table[t][leaf_name] = None
            
    return pd.DataFrame.from_dict(ancestor_table, orient='index')

 
def initialize_llh_objects(
    genetrees: List[dendropy.Tree], tip_traits: Dict, max_time: float
):
    llhobjs = [tree_llh_continuous(gt) for gt in genetrees]
    # compute tree heights and initialize pruning
    heights = [max(obj.tree.calc_node_root_distances()) for obj in llhobjs]
    horizon = max(max_time, max(heights))
    for obj in llhobjs:
        obj.initialize_pruning(tip_traits, horizon)
    return llhobjs


def compute_total_llh(llhobjs, rate_vec: Iterable[float]) -> float:
    # sum over trees, respecting individual weights if your calc_llh_per_tree uses them internally
    return float(sum(calc_llh_per_tree((rate_vec, obj)) for obj in llhobjs))


def compute_ancestral_state_dfs(llhobjs, time_grid: np.ndarray):
    return [ancestral_state_map(obj.tree, time_grid) for obj in llhobjs]


def build_state_matrix(
    dfs: List[pd.DataFrame],
    tree: dendropy.Tree,
    time_grid: np.ndarray,
    exclude_leaves: Optional[List[str]] = None,
    n_bins: int = 100,
) -> np.ndarray:
    exclude = set(exclude_leaves or [])
    leaves = {
        leaf.taxon.label: leaf
        for leaf in tree.leaf_node_iter()
        if leaf.taxon.label not in exclude
    }
    # rows: time (ascending); cols: trait bins (fixed 100)
    mat = np.zeros((len(time_grid), n_bins), dtype=float)

    for t_idx, t in enumerate(time_grid):
        # collect state distributions for each kept leaf at time t across all dfs
        row_accum = []
        for name in leaves.keys():
            # some dfs might miss exact t if slicing/coercion fails; skip safely
            states = []
            for df in dfs:
                try:
                    s = df.loc[t, name].state  # expected to be a 1D probability vector
                    states.append(s)
                except Exception:
                    continue
            if not states:
                continue
            # average normalized states across gene trees
            states = [np.array(s, dtype=float) for s in states]
            # normalize each, guard zero
            normed = []
            for s in states:
                ssum = s.sum()
                normed.append(s / ssum if ssum > 0 else s)
            mean_state = np.mean(normed, axis=0)
            row_accum.append(mean_state)

        if row_accum:
            mat[t_idx, :] = np.sum(row_accum, axis=0)

    # reverse in time so earlier at bottom
    return mat[::-1, :]



# Plotting
def plot_state_contour(
    state_matrix: np.ndarray,
    trait_vec: np.ndarray,
    time_grid: np.ndarray,
    lower: int,
    higher: int,
    ax: Optional[plt.Axes] = None,
    cmap: str = "plasma",
    fontsize: int = 20,
    colorbar: bool = False,
) -> plt.Axes:
    sns.set_theme(style="white")
    ax = ax or plt.gca()

    # slice columns and arrange coordinates
    state = state_matrix[:, lower:higher]
    x = np.asarray(trait_vec[lower:higher])
    y = np.asarray(time_grid[::-1])
    X, Y = np.meshgrid(x, y)

    vmin, vmax = 0.0, 0.3
    levels = np.linspace(vmin, vmax, 50)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    cs = ax.contourf(X, Y, state, levels=levels, cmap=cmap, norm=norm)

    ax.set_xlabel("Trait", fontsize=fontsize, labelpad=20)
    ax.set_ylabel("Time", fontsize=fontsize)
    ax.tick_params(axis="x", labelsize=fontsize)
    ax.tick_params(axis="y", labelsize=fontsize)

    if colorbar:
        cbar = plt.colorbar(cs, ax=ax)
        cbar.ax.tick_params(labelsize=fontsize - 2)
        cbar.locator = ticker.MaxNLocator(nbins=5)
        cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        cbar.update_ticks()

    return cs


def ancestral_state_pipeline(
    cfg: Dict, 
    prefix: str = "example",
    treetype: str = "st",
    trait_ids: Optional[List[int]] = None,
    lower: int = 0,
    higher: int = 100,
    ax: Optional[plt.Axes] = None,
    add_colorbar: bool = False,
):
    final_time = cfg["final_time"]
    num_slicing = cfg["num_slicing"]
    speciation_time = cfg["speciation_time"]
    trait_ids = trait_ids or cfg.get("traits", [])
    base_dir = cfg["base_dir"]
    time_grid = np.linspace(0, final_time, num_slicing)

    paths = build_paths(prefix, treetype, base_dir)
    sptree, ratetree, segment2rate = load_species_and_rate_trees(paths)
    genetrees = load_gene_trees(paths, sptree.taxon_namespace)

    # prep gene trees
    slice_and_map_gene_trees(sptree, genetrees, segment2rate, time_grid)

    # load traits and estimated rate 
    tip_traits_list, est_rate = load_traits_and_rates(paths)

    # one figure per trait_id
    for trait_idx in trait_ids:
        tip_trait = tip_traits_list[trait_idx]
        llhobjs = initialize_llh_objects(genetrees, tip_trait, max_time=time_grid.max())
        curr_rate = est_rate.loc[est_rate.trait_idx == trait_idx, "rate_0"].values
        _total_llh = compute_total_llh(llhobjs, curr_rate)

        # ancestral states per tree
        dfs = compute_ancestral_state_dfs(llhobjs, time_grid)

        # pull representative tree and its trait vector
        ref_tree = llhobjs[0].tree
        trait_vec = llhobjs[0].trait_vec

        # build matrix (sum over leaves except excluded)
        # 0-T2
        curr_to_t2 = build_state_matrix(
            dfs=dfs,
            tree=ref_tree,
            time_grid=time_grid[time_grid < speciation_time],
            exclude_leaves=["A2", "A3"],
            n_bins=100,
        )

        # T2-deep back
        t2_back = build_state_matrix(
            dfs=dfs,
            tree=ref_tree,
            time_grid=time_grid[time_grid >= speciation_time],
            exclude_leaves=["A3"],
            n_bins=100,
        )

        matrix = np.concat((t2_back, curr_to_t2))

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6), dpi=100)
        cs = plot_state_contour(
            matrix,
            trait_vec,
            time_grid,
            lower,
            higher,
            ax=ax,
            cmap="plasma",
            colorbar=add_colorbar,
        )

        return ax, cs
