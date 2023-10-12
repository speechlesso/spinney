from utils import *
import argparse, os, time, re
from multiprocessing import Pool
import threading
file_lock = threading.Lock()

NUM_PROCESS = 5

def opt_func_pool(rates, llhobjs):
    with Pool(NUM_PROCESS) as pool:
        llhs = pool.map(calc_llh_per_tree, [(rates, x) for x in llhobjs])
    return sum(llhs)

def opt_func(rates, llhobjs):
    llhs = [ calc_llh_per_tree((rates,x)) for x in llhobjs]
    return sum(llhs)


def process_item(tip_trait, trees, outfile, num, condition):
    init_time = time.time()

    # init
    llhobjs = [tree_llh_continous(x) for x in trees]
    trees_height = [max(llhobj.tree.calc_node_root_distances()) for llhobj in llhobjs]
    cov_mat = sum([llhobj.C_mat() * llhobj.tree.weight for llhobj in llhobjs])
    if llhobjs[0].tree.weight == 1:
        # if input trees are sampled gene trees
        cov_mat /= len(llhobjs)

    [x.initialize_pruning(tip_trait, max(trees_height)) for x in llhobjs]
    init_val = np.var([v for v in tip_trait.values()]) / np.max(trees_height)

    # number of rate
    numrate = len( np.unique([v for v in segment2rate.values()]) )
    init_rates = (init_val,) * numrate

    # opt
    mle = mle_sigma2(sorted(trees[0].leaf_nodes(), key=lambda x: x.taxon.label), tip_trait, cov_mat)
    res = minimize(opt_func_pool,
                   init_rates,
                   args=(llhobjs),
                   method='Nelder-Mead',
                   bounds=((1e-8, np.inf),) * numrate,
                   #constraints=[{"type": "ineq", "fun": opt_func}],
                   #options={}
                   )
    res2file = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
    len(trees), res.success + 0, res.fun, '\t'.join(['%s'%x for x in res.x]), mle,\
    tip_trait[llhobjs[0].leaves[0].taxon.label], condition, num,\
    round(time.time() - init_time))
    print(res2file)
    with file_lock:
        with open(outfile, 'a') as file:
            file.write(res2file + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='traits saved as csv', type=str)
    parser.add_argument('-N', help='number of tree for optimization', type=int, default=10)
    parser.add_argument('-genetrees', help='file for gene tree of optimization. branch in 2N unit', type=str, default='')
    parser.add_argument('-speciestree', help='file for species tree. branch in 2N unit', type=str, default='')
    parser.add_argument('-ratetree', help='rate family file', type=str, default=None)
    parser.add_argument('-condition', help='comments',type=str,default='')
    parser.add_argument('-timeslice', help='add internal unifurcated nodes', type=int, default=None)

    args = parser.parse_args()

    ## DEBUG
    # args.f = './trait/4taxa_high.csv'
    # args.N = 10
    # args.genetrees = './tree/4taxa_high.trees.tmp'
    # args.speciestree = './tree/4taxa_high.sptree.nwk'
    # args.ratetree = './tree/4taxa_high.rtree.nwk'
    # args.condition = '4taxa_1p_gt'

    # TREE
    nwkstring = read_nwk_trees(args.genetrees)
    sptree = nwk2tree( read_nwk_trees(args.speciestree)[0] )
    # add time slicing (speciation time) node for ancestral state
    speciation_time = [node.distance_from_tip() for node in sptree.nodes() if len(node.child_nodes()) == 2]
    speciation_time = list( set(speciation_time) )
    speciation_time.sort(reverse=True)
    # exclude root node
    # speciation_time = speciation_time[:-1]
    # speciation_time = []

    ratetree = dendropy.Tree.get(data=read_nwk_trees(args.ratetree)[0],
                                 schema="newick",
                                 taxon_namespace=sptree.taxon_namespace)
    [setattr(edge, 'sp_segment', idx) for idx, edge in enumerate(sptree.preorder_edge_iter())]
    segment2rate = {sedge.sp_segment: int(redge.length) \
                    for sedge, redge in \
                    zip([se for se in sptree.preorder_edge_iter()], [re for re in ratetree.preorder_edge_iter()])}

    # read newick gene trees with weight
    tree_weights = parse_trees(args.genetrees)
    genetrees = []
    for t, freq in tree_weights.items():
        genetrees.append( dendropy.Tree.get(data=t, schema="newick", taxon_namespace=sptree.taxon_namespace) )
        genetrees[-1].weight = freq

    args.N = len(genetrees)

    # TRAIT
    sim_ancester = pd.read_csv(args.f)
    traits = sim_ancester
    tip_traits = [traits.loc[num].to_dict() for num in range(traits.shape[0])]


    outfile = './result/%s.tsv'%args.condition

    expect_trait = lambda x, y: sum(x.trait_vec * y.state) / sum(y.state)
    for num in range(50):
        # np.random.seed(num)
        # trees = [ nwk2tree(x) for x in np.random.choice(nwkstring, size=args.N)]
        # trees = [dendropy.Tree.get(data=newick_tree, schema="newick", taxon_namespace=sptree.taxon_namespace) \
                  # for newick_tree in np.random.choice(nwkstring, size=args.N)]
        for treeno, t in enumerate(genetrees):
            [time_slice_node( t, 1.001 * timeslice ) for timeslice in speciation_time]
            segment_mapping(sptree, t, speciation_time)
            for node in t.nodes():
                if hasattr(node.edge, 'sp_segment'):
                    setattr(node, 'rate_family', segment2rate[node.edge.sp_segment]-1)
                # node.comments = segment2rate[node.edge.sp_segment] - 1
            # t.weight = 1/len(trees)

            for idx, x in enumerate(t.preorder_node_iter()):
                try:
                    x.label += '-rt_%s'%(str( x.rate_family ))
                except AttributeError:
                    x.label += '-None'

            #print(treeno, t.as_string('newick'))

        # INPUT SINGLE TRAIT
        process_item(tip_traits[num], trees, outfile, num, args.condition)
