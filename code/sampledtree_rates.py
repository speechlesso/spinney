# load packages
from utils import *
import argparse, os, time, re
from multiprocessing import Pool

def opt_func_pool(rates, llhobjs):
    with Pool(NUM_PROCESS) as pool:
        llhs = pool.map(calc_llh_per_tree, [(rates, x) for x in llhobjs])
    return sum(llhs)


def opt_func(rates, llhobjs):
    llhs = [ calc_llh_per_tree((rates,x)) for x in llhobjs]
    return sum(llhs)


def rate_inference(tip_trait, trees, outfile, num, condition):
    init_time = time.time()

    # initialization
    llhobjs = [tree_llh_continuous(x) for x in trees]
    tree_heights = [max(llhobj.tree.calc_node_root_distances()) for llhobj in llhobjs]
    [x.initialize_pruning(tip_trait, max(tree_heights)) for x in llhobjs]

    # rates
    numrate = len( np.unique([v for v in segment2rate.values()]) )
    init_value = np.var([v for v in tip_trait.values()]) / np.max(tree_heights)
    init_rates = (init_value,) * numrate

    # inference with c star
    cstar_mat = sum([llhobj.C_mat() * llhobj.tree.weight for llhobj in llhobjs])
    cstar_est = mle_sigma2(sorted(trees[0].leaf_nodes(), key=lambda x: x.taxon.label), tip_trait, cstar_mat)

    # inference
    res = minimize(opt_func_pool,
                   init_rates,
                   args=(llhobjs),
                   method='Nelder-Mead',
                   bounds=((1e-8, np.inf),) * numrate,
                   )
    
    outputs = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
    len(trees), res.success + 0, res.fun, \
    '\t'.join(['%s'%x for x in res.x]), cstar_est,\
    tip_trait[llhobjs[0].leaves[0].taxon.label], condition, num,\
    round(time.time() - init_time))

    print(outputs)
    with open(outfile, 'a') as file:
        file.write(outputs + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='traits saved as csv', type=str)
    parser.add_argument('-N', help='number of tree for optimization', type=int, default=10)
    parser.add_argument('-genetrees', help='file for gene tree of optimization. branch in 2N unit', type=str, default='')
    parser.add_argument('-speciestree', help='file for species tree. branch in 2N unit', type=str, default='')
    parser.add_argument('-ratetree', help='rate family file', type=str, default=None)
    parser.add_argument('-condition', help='comments',type=str,default='')
    parser.add_argument('-timeslice', help='add internal unifurcated nodes', type=int, default=None)
    parser.add_argument('-num_process', help='Number of process for parallelization', type=int, default=5)

    args = parser.parse_args()
    global NUM_PROCESS
    NUM_PROCESS = args.num_process

    # TREE
    nwkstring = read_nwk_trees(args.genetrees)
    sptree = nwk2tree( read_nwk_trees(args.speciestree)[0] )
    # add  speciation time (time slice) node for ancestral state
    speciation_time = [node.distance_from_tip() for node in sptree.nodes() if len(node.child_nodes()) == 2]
    speciation_time = list( set(speciation_time) )
    speciation_time.sort(reverse=True)
    # exclude root node
    # speciation_time = speciation_time[:-1]
    # speciation_time = []

    # construct a tree specifying rate
    # rate tree has identical topology as species tree
    ratetree = dendropy.Tree.get(data=read_nwk_trees(args.ratetree)[0],
                                 schema="newick",
                                 taxon_namespace=sptree.taxon_namespace)
    [setattr(edge, 'sp_segment', idx) for idx, edge in enumerate(sptree.preorder_edge_iter())]
    segment2rate = {sedge.sp_segment: int(redge.length) \
                    for sedge, redge in \
                    zip([se for se in sptree.preorder_edge_iter()], [re for re in ratetree.preorder_edge_iter()])}
    
    # read gene trees with weight
    tree_weights = parse_trees(args.genetrees)
    genetrees = []
    for nwk, freq in tree_weights.items():
        genetrees.append( dendropy.Tree.get(data=nwk, 
                                            schema="newick", 
                                            taxon_namespace=sptree.taxon_namespace) )
        genetrees[-1].weight = freq
    
    args.N = len(genetrees)
    # Assign rate to each segment in all gene trees
    for _, t in enumerate(genetrees):
        [time_slice_node( t, 1.00001 * timeslice ) for timeslice in speciation_time]
        segment_mapping(sptree, t, speciation_time)
        for node in t.nodes():
            if hasattr(node.edge, 'sp_segment'):
                setattr(node, 'rate_family', segment2rate[node.edge.sp_segment]-1)

    # TRAIT
    traits = pd.read_csv(args.f, comment = "#", dtype=float)
    tip_traits = [traits.loc[num].to_dict() for num in range(traits.shape[0])]


    # INFERENCE
    outfile = './%s.tsv'%args.condition
    for trait_idx in range(100):
        rate_inference(tip_traits[trait_idx], genetrees, outfile, trait_idx, args.condition)
