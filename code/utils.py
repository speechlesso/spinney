import numpy as np
import pandas as pd
from Bio import Phylo
from io import StringIO
from scipy.optimize import minimize
from itertools import product, combinations
import math
from math import pi
import copy, pickle, re
import dendropy, ete3


#%% TREE process
def read_nwk_trees(treefile):
    # read newick trees from a file
    nwks = []
    with open (treefile,'r') as f:
        for line in f.readlines():
            if line!='\n' and len(re.findall("\(", line))!=0:
                nwks.append(line.split('\n')[0])
    return nwks

def nwk2tree(newick_tree, is_4N_unit = False):
    tree = dendropy.Tree.get(data=newick_tree,schema="newick")
    #numInternalNode = 1
    internal_node = [x for x in tree.postorder_internal_node_iter()]
    for clade in internal_node:
        #clade.name = 'N' + str(numInternalNode)
        #numInternalNode += 1
        if clade.edge_length is None:
            clade.edge_length = 0.0
    if is_4N_unit:
        for x in tree.postorder_node_iter():
            x.edge_length *= 2
    return tree

# Add rate as attribute for each branch segment with a rate tree
def branch_rate_family(treeobj, sptree, rate_tree=None):
    if rate_tree != None: # multiple rate
    # rate_tree has identical topology as speceis tere
        # for c1, c2 in zip(treeobj.nodes(), rate_tree.nodes()):
            # setattr(c1, 'rate_family', int(c2.edge.length)-1)
        [setattr(sedge, 'rate_family', int(redge.length)) \
         for sedge, redge in
         zip([se for se in sptree.preorder_edge_iter()], [re for re in rate_tree.preorder_edge_iter()])]
        # TODO: gene tree segment ->segment_mapping id -> rate_family index
    else:
    # single rate
        [setattr(clade, 'rate_family', 0) for clade in treeobj.nodes()]

    # [setattr(edge, 'sp_segment', idx) for idx, edge in enumerate(self.sptree.preorder_edge_iter())]
    # [setattr(sedge, 'rate', redge.length) \
    #  for sedge, redge in
    #  zip([se for se in self.sptree.preorder_edge_iter()], [re for re in self.rate_tree.preorder_edge_iter()])]
    return treeobj

def parse_trees(treefile):
    with open(treefile, 'r') as f:
        lines = f.readlines()
    # Remove any whitespace or newline characters from each line
    lines = [line.strip() for line in lines]

    # Check where the gene trees and genetree freqs sections start and end
    begin_genetrees = lines.index("begin genetrees;")
    end_genetrees = lines.index("end;")
    begin_genetreefreqs = lines.index("begin genetreefreqs;")

    # Extract gene trees and their frequencies
    genetreesnwk = lines[begin_genetrees + 1:end_genetrees]
    genetreefreqs = lines[begin_genetreefreqs + 1:-1]  # -1 to exclude the "end;" line

    # Convert frequencies from string to float
    genetreefreqs = [float(freq) for freq in genetreefreqs]

    tree_weights = dict(zip(genetreesnwk, genetreefreqs))

    return tree_weights

def get_speciation_time(tree):
    brs = [x.edge_length for x in tree.leaf_nodes()]
    T1 = min(brs)
    T2 = max(brs)
    return T1, T2

def observed_discordance(mstree):
    if type(mstree) is str:
        nwks = read_nwk_trees(mstree)
    else:
        nwks = mstree

    if nwks[0].startswith('[&R]'):
        nwks = [x.split('[&R] ')[1] for x in nwks]

    topos = [extract_topology(x) for x in nwks]
    topos = [extract_topology(x) for x in nwks]
    topo, cnt = np.unique(topos, return_counts=True)
    print('The majority topology is %s %.2f%%'%(topo[cnt.argmax()], 100*max(cnt)/len(topos)))

# define Python user-defined exceptions
class IncorrectTreeInput(Exception):
    pass

# parse topology in newick tree
def extract_topology(newick):
    topology = re.sub(':[0-9]*\.[0-9]*','',newick)
    return topology

# Add speciation time slice as internal nodes to gene trees
def time_slice_node(tree, slice2tip):
    root2tip = tree.max_distance_from_root()
    root2slice = root2tip - slice2tip
    if root2slice < 0 or slice2tip < 0:
        print('distance for root to time slice is %.2f, distance from time slice to tip is %.2f' % (
        root2slice, slice2tip))
        return tree
        #raise BaseException('Time slice is greater than the height of tree.')

    # get possible parents for time-slice node in each lineage
    parent4slice = [x for x in tree.levelorder_node_iter() if x.distance_from_tip() >= slice2tip]

    # add time slice as new nodes
    for parent_node in parent4slice:
        slice2parent = parent_node.distance_from_tip() - slice2tip
        children = parent_node.child_nodes()
        if len(children) > 1 :
            for child in children:
                new_edge = child.edge_length - slice2parent
                if new_edge >= 0:
                    child = parent_node.remove_child(child)
                    child.edge_length = new_edge
                    slice_node = parent_node.new_child(edge_length=slice2parent)
                    slice_node.add_child(child)
                    setattr(slice_node, 'time_slice', True)
                    tree.update_bipartitions(suppress_unifurcations=False)

    for idx,x in enumerate(tree.preorder_internal_node_iter()):
        x.label = 'N%s'%idx
        try:
            x.time_slice
        except AttributeError:
            setattr(x, 'time_slice', False)

    for idx,x in enumerate(tree.leaf_node_iter()):
        x.label = x.taxon.label
        setattr(x, 'time_slice', True)

    # [setattr(x, 'time_slice', True) for x in tree.leaf_nodes()]
    # print(tree.as_string('newick'))
    return tree

def find_path_to_root(target):
    path = [target.edge]
    node = target
    while node.parent_node is not None:
        node = node.parent_node
        path.append(node.edge)
    return path

def rename_internal(tree):
    numInternalNode = 1
    internal_node = [x for x in tree.postorder_internal_node_iter()]
    internal_node.reverse()
    for clade in internal_node:
        clade.label = 'N' + str(numInternalNode)
        numInternalNode += 1
    return tree

def rename_leaf(tree):
    for node in tree.leaf_nodes():
        node.label = node.taxon.label

def get_3taxa_newick(N, t1=4000, t2=50000):
    br1 = t1 / (2 * N)
    br2 = (t2 - t1) / (2 * N)
    newick_string = '[&R] ((sp1:%s,sp2:%s):%s, sp3:%s);'%(br1, br1, br2, br1 + br2)
    return newick_string


def segment_mapping(sptree, genetree, speciation_time):
    '''
    Assign rate from the species tree to each segment of a gene tree
    '''
    # speciation_time.sort()
    # input species tree, map each edge to an integer
    leaf_names = [x.taxon.label for x in sptree.leaf_nodes()]
    leaf_names.sort()
    sedges = [edge for edge in sptree.preorder_edge_iter()]
    sedge2int = {edge:idx for idx, edge in enumerate(sedges)}
    # sedge2interval = {edge:(edge.head_node.distance_from_tip(), edge.tail_node.distance_from_tip(), \
    #                             [x.taxon.label for x in edge.head_node.leaf_nodes()]) \
    #                   for edge in sedges[:-1]}

    # get path from root to edge, each element is an edge object in dendropy
    gpaths = [ find_path_to_root(leaf) for leaf in genetree.leaf_node_iter() ]

    for gpath in gpaths:
        leaf_taxon = gpath[0].head_node.taxon.label
        spath = [ find_path_to_root(leaf) \
                  for leaf in sptree.leaf_node_iter()\
                  if leaf.taxon.label == leaf_taxon][0]
        # edge x : root------ tail ---x--- head
        for edge in gpath:
            try:
                # avoid redundant assigning
                edge.sp_segment
            except AttributeError:
                # find mapped branch in sptree
                idx = [edge.head_node.distance_from_tip() >= t for t in speciation_time]

                if any(idx):# non-terminal edge
                    distance2tip = max([ t for t, flag in zip(speciation_time, idx) if flag ])
                else:  # terminal edge
                    distance2tip = 0.0
                
                mapped_sedge = [se for se in spath if se.head_node.distance_from_tip()  <= distance2tip]
                setattr( edge, 'sp_segment', sedge2int[mapped_sedge[-1]] )
    # root node
    # [edge for edge in genetree.preorder_edge_iter()][0].sp_segment = -1
    return genetree


# TODO: NOT WORKING
def topology_test(tree, query_tree, symmetric_distance=False):
    sub_tree_bitmask = tree.taxon_set.get_taxa_bitmask(labels=query_tree.taxon_set.labels())
    sub_tree_mrca_node = tree.mrca(leafset_bitmask=sub_tree_bitmask)
    sub_tree_newick = sub_tree_mrca_node.as_newick_string()
    sub_tree_ = dendropy.Tree()
    sub_tree_.read_from_string(sub_tree_newick, schema='newick')
    sd = sub_tree_.symmetric_difference(query_tree)
    if symmetric_distance == True:
        return sd

    elif sd == 0:
        return True

    else:
        return False

def common_ancestor(tree, targets):
    """Most recent common ancestor (clade) of all the given targets.
       Edge cases:
        - If no target is given, returns self.root
        - If 1 target is given, returns the target
        - If any target is not found in this tree, raises a ValueError
    Modified from Bio.Phylo.BaseTree.common_ancestor for dendropy tree object
    """
    paths = [[ x.head_node for x in find_path_to_root(t)[:-1]] for t in targets]
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


def is_unifurcation(node):
    return len(node.child_nodes()) == 2


#%% three-taxa tree constructor
# construct gene trees for 3-taxa species tree
class build_gene_tree():
    def __init__(self, nwk_st, is_4N_unit = False):
        self.is_4N_unit = is_4N_unit
        self.tree = nwk2tree(nwk_st, is_4N_unit)

    # works for a three taxa case
    def construct_gene_tree(self, add_timeslice=False):
        tipname = [] # [tip2inner, tip2inner, tip2root ]
        # read branch length in species tree
        # speices tree: ((A:x, B:x)N1:z, C:y);
        # get speciation time in species tree
        T1, T2 = get_speciation_time(self.tree)

        for node in self.tree.leaf_nodes():
            parent = node.parent_node
            if parent.distance_from_root() != 0:
                # leave node connected to internal node
                x = node.edge.length
                tipname.append(node.taxon.label)
            else:
                # leave node connected to root
                y = node.edge.length
                tipname.append(node.taxon.label)

        # internal node connected to root
        z = y - x

        if y != x+z:
            # print('Incorrect tree input!')
            raise IncorrectTreeInput

        # LS gene tree 1: ((A:m, B:n)l, C:k);
        m = x + ( 1 - z/(np.exp(z)-1) )
        n = x + ( 1 - z/(np.exp(z)-1) )
        k = x + z + 1
        l = k - m
        n1, n2, n3 = tipname
        gt1_nwk = '[&R] ((%s:%.10f, %s:%.10f):%.10f, %s:%.10f);' %(n1, m, n2, n, l, n3, k)
       
        # calc repeated values
        k = tip2root = (x + z) + 4/3
        m = tip2inter= (x + z) + 1/3
        inter = 1.0

        # ILS trees
        gt2_nwk = '[&R] ((%s:%.10f, %s:%.10f):%.10f, %s:%.10f);' %(n1, tip2inter, n2, tip2inter, inter, n3, tip2root)
        gt3_nwk = '[&R] ((%s:%.10f, %s:%.10f):%.10f, %s:%.10f);' %(n2, tip2inter, n3, tip2inter, inter, n1, tip2root)
        gt4_nwk = '[&R] ((%s:%.10f, %s:%.10f):%.10f, %s:%.10f);' %(n1, tip2inter, n3, tip2inter, inter, n2, tip2root)

        gts = []

        # weigths
        concordant_freq = 1-np.exp(-z)
        discordant_freq = np.exp(-z)/3


        weights = [concordant_freq,] + [discordant_freq,]*3
        for gt_nwk, w, ids in zip([gt1_nwk, gt2_nwk, gt3_nwk, gt4_nwk], weights, ['ls','ils','ils','ils']):
            gt = nwk2tree(gt_nwk)
            gt.weight = w
            gt.label = ids
            if add_timeslice:
                gt = time_slice_node(gt,T1)
                gt = time_slice_node(gt,T2)
            gts.append( gt )
        return gts
    

#%% Class for calculating likelihood for continous trait
class tree_llh_continuous():
    def __init__(self, treeobj, vec_size=100):
        self.tree = copy.deepcopy(treeobj)
        self.tree.root = self.tree.nodes()[0]
        self.leaves = sorted(self.tree.leaf_nodes(), key= lambda x: x.taxon.label)
        self.ordered_vertices = [x for x in self.tree.postorder_node_iter()]
        
        self.node_to_num = {n:i for i, n in \
                            enumerate([x for x in self.tree.postorder_node_iter() ])}
        self.num_to_node = {v : k for k, v in self.node_to_num.items()}
        self.stationary = np.ones(vec_size)/vec_size

    def C_mat(self):
        # construct covariance matrix
        # as non censored method in O'Meara 2006 

        # diag elements
        leaves = sorted(self.tree.leaf_nodes(), key=lambda x: x.taxon.label)
        cov = np.diag( self.tree.calc_node_root_distances() )
        # off-diag elements
        tip_to_num = {n.taxon:i for i, n in enumerate(leaves)}
        leaves_pair = [(t1, t2) for t1,t2 in combinations(leaves, 2)]

        for pair in leaves_pair:
            taxon1, taxon2 = pair
            mrca = common_ancestor(self.tree, pair)
            r, c = [tip_to_num[x.taxon] for x in pair]
            cov[r,c] = cov[c,r] = mrca.distance_from_root()
        return cov


    def initialize_pruning(self, tip_trait, tree_height, vec_size = 100):
        '''
        INITIALIZE FOR ONE TRAIT EVERY TIME!
        THE TRAIT RANGE VARIES DEPENDING ON TIP TRAITS!

        tip_trait:a dict {leavename:single_trait_value}
        vec_size: vector length for range of trait
        '''

        # discretize continous traits
        max_val = max(  ([v for v in tip_trait.values()])  )
        min_val = min(  ([v for v in tip_trait.values()])  )
        dev =  np.sqrt( tree_height )
        self.trait_vec = np.linspace(min_val - 10*dev, max_val + 10*dev, vec_size).flatten()

        # calculate state of tips
        tip_state = {}
        for key, val in tip_trait.items():
            temp_vec = np.zeros(self.trait_vec.size)
            # trait value matches one in discretized vector exactly
            val_index = np.argwhere(self.trait_vec==val).flatten()
            if len(val_index) != 0:
                temp_vec[val_index] = 1.0
            else:
                # trait value falls btween two numbers corresponding to trait range, e.g. 0.5 in [0,1]
                upperind = np.min( np.argwhere(self.trait_vec>val) )
                lowerind = upperind -1
                stepsize = self.trait_vec[upperind] - self.trait_vec[lowerind]
                weight = (val - self.trait_vec[lowerind])/stepsize
                temp_vec[upperind] = weight
                temp_vec[upperind-1] = 1 - temp_vec[upperind]
            tip_state.update({key:temp_vec})

        for leaf in self.tree.leaf_nodes():
            setattr(leaf, 'state', tip_state[leaf.taxon.label])


    def calc_llh_with_tree(self, rate):
        nstate = self.stationary.size
        # initialize state for internal nodes
        for clade in self.tree.internal_nodes():
            setattr(clade, 'state', np.zeros(nstate))

        # cache probability density matrix for trees
        for clade in self.ordered_vertices:
            if clade == self.tree.root:
                setattr(clade, 'probmat', np.eye(nstate))
            else:
                prob_mat = np.zeros((nstate, nstate))
                for row in range(nstate):
                    for col in range(row+1):
                        prob_mat[row,col] = bm_prob_density(self.trait_vec[col], 
                                                            self.trait_vec[row], 
                                                            rate[clade.rate_family], 
                                                            clade.edge_length)
                        prob_mat[col,row] = prob_mat[row,col]
                row_sum = prob_mat.sum(axis = 1)
                setattr(clade, 'probmat', prob_mat/row_sum)

        # pruning algorithm
        for parent in self.tree.postorder_internal_node_iter():
            child_llh = [child.probmat @ child.state for child in parent.child_nodes() ]
            if len(child_llh) == 1:
                parent.state = child_llh[0]
            else:
                # consider two children at most
                parent.state = child_llh[0] * child_llh[1]

        negative_llh =  - np.log( self.stationary.dot(self.tree.root.state)  )
        return negative_llh


    def expected_trait(self):
        for clade in self.tree.nodes():
            trait_temp = sum((self.trait_vec * clade.state))/sum(clade.state)
            setattr(clade, 'trait', trait_temp)
            print(clade, trait_temp)



#%% HELPER FUNCTION FOR LIKLIHOOD
def bm_prob_density(x, x0, sigma2, t):
    # x:trait value; x0: initial value; sigma2: rate; t:branch length
    t = t if t!=0 else 1e-8 # avoid divided by zero
    coef1 = 1/((2*pi*t*sigma2)**0.5)
    coef2 = -(x-x0)**2/(2*sigma2*t)
    prob = coef1*np.exp(coef2)
    return prob


def mle_sigma2(leaves, tip_trait, cov):
    leaves_name = [x.taxon.label for x in leaves]
    n_taxa = len(tip_trait)
    ones = np.ones([n_taxa,1])
    x = np.array([tip_trait[x] for x in leaves_name]).reshape(n_taxa,1)
    c = cov
    inv_c = np.linalg.inv(c)
    a = (np.linalg.inv(ones.T @ inv_c @ ones)) @ (ones.T @ inv_c @ x)
    #a = ones @ (ones.T @ inv_c @ ones) @ (ones.T @ inv_c @ x)
    return ( (x-a).T @ inv_c @ (x-a)/n_taxa )[0][0]


def llh_surface(tip_trait, cov, sigma2):
    n_taxa = len(tip_trait)
    sort_key = [k for k in tip_trait.keys() ]
    sort_key.sort()
    x = np.array([tip_trait[k] for k in sort_key]).reshape(n_taxa,1)

    inv_c = np.linalg.inv(cov * sigma2)
    det_c = np.linalg.det(cov * sigma2)

    denom = np.exp(-0.5 * (x - np.mean(x)).T @ inv_c @ (x - np.mean(x)))
    nome = np.sqrt((2*np.pi)**n_taxa * det_c)
    return -np.log(denom/nome).flatten()[0]


def calc_llh_per_tree(x):
    rates, llhobj  = x
    return llhobj.calc_llh_with_tree(rates) * llhobj.tree.weight


def init_llh_per_tree(x):
    llhobj, tip_trait, tree_height = x
    llhobj.initialize_pruning(tip_trait, tree_height)
    return llhobj


#%% continous trait simulator
def get_continous_trait(mean, cov_mat, seed=None):
    # construct cov_var matrix then sample trait for the tips.
    if seed != None:
        np.random.seed(seed)
    return np.random.multivariate_normal(mean, cov_mat)


# simulate continous trait with brownian motion
def sim_trait_val(species_tree, rate, root_trait, seed=None):
    # TODO: multiple rates for construct cov
    leaves = species_tree.get_terminals()
    n_taxa = len(leaves)
    # diag elements
    diag = species_tree.distance(leaves[0])
    cov = np.diag( [diag,]*n_taxa )

    # off-diag elements
    internal = [x for x in species_tree.get_nonterminals() if x!=species_tree.root ]
    tip_to_num = {n.name:i for i, n in enumerate(leaves)}

    leaves_pair = [x for x in combinations(leaves, 2)]
    for pair in leaves_pair:
        ca = species_tree.common_ancestor(pair)
        r, c = [ tip_to_num[x.name] for x in pair ]
        cov[r,c] = cov[c,r] = rate * species_tree.distance(ca)
    tip_init = dict(  zip([x.name for x in leaves], get_continous_trait(root_trait, cov, seed))  )  #{leafname:trai_value}
    #print('cov4simulating trait\n', cov)
    return tip_init

def sim_ancestral_state(tree, rate, root_dist=None, seed=None):
    np.random.seed(seed)
    if root_dist == None:
        root_dist = 0

    for parent in tree.find_clades():
        if parent == tree.root:
            setattr(parent, 'simtrait', root_dist)
        for child in parent.clades:
            trait_temp = np.random.normal(loc = parent.simtrait, 
                                          scale = np.sqrt(child.branch_length * rate[child.rate_family]) )
            setattr(child, 'simtrait', trait_temp)
'''
#### TEST FOR sim_ancestral_state with single rate
t = nwk2tree("[&R] ((sp1:2,sp2:2):2,sp3:4);")
time_slice_node(t, 3)
assign_branch_rate(t, 0.5)
sim_ancestral_state(t,root_dist = 10)
'''


#%% OTHER HELPER FUNCTION
def g(i,j,T):
    prob = 0
    for k in range(j, i+1):
        upper =  math.factorial(j+k-2) *  math.factorial(i) / ( math.factorial(j-1) * math.factorial(i-k))
        lower =  math.factorial(j) * math.factorial(k-j) *  math.factorial(i+k-1) /  math.factorial(i-1)
        prob += np.exp(-k*(k-1)*T/2) * (2*k-1)*(-1)**(k-j) * upper/lower
    return prob

def four_taxa_concordance(T2, T3):
    # Species Tree Topology Is (((A:t4,B:t4):t3,C:t3+t4):t2,D:t2+t3+t4);
    # function is derived from Noah A. Rosenberg 2002 TABLE V
    prob = g(2,1,T3)*( g(2,1,T2) + g(2,2,T2)/3) \
           + g(2,2,T3)*(g(3,1,T2)/3 + g(3,2,T2)/9 + g(3,3,T2)/18)
    return prob

def branch_scale_coef(discord=1e-3):
    N1 = effect_N(1e-3)
    N2 = effect_N(discord)
    return N2/N1

def effect_N(d=1e-3):
    # br = internal branch; d = level of discordance
    if d>2/3:
        raise BaseException('The discordance should be less than 2/3')
    return -0.5/(np.log(3*d/2))

# threee-taxa case for different discordance
# adjust N
def sptree_4taxa(N=2000, t4=4000, t3=50000, t2=90000):
    x = t4/(2*N)
    y = t3/(2*N)
    z = t2/(2*N)
    disc = 1 - four_taxa_concordance(z-y, y-x)
    #print('The discordance is {:.3f}% with N = {:0d}'.format(dis*100, N))
    newick_tree = '(((sp1:%.4f,sp2:%.4f):%.4f,sp3:%.4f):%.4f,sp4:%.4f);'%(x,x,y-x,y,z-y,z)
    return newick_tree, disc

def calc_stat(state, traitval):
    mean = np.sum(state*traitval)
    mode = traitval[np.argmax(state)]
    var = np.sum(state*traitval**2) - mean**2
    return mean, mode, var


#%% DEBUG

# sptreenwk = '[&R] ((sp2:0.02500,sp3:0.02500):0.28750,sp4:0.31250);'
# sptree = nwk2tree(sptreenwk)
# gen_gt = build_gene_tree(sptreenwk)
# trees_ts = gen_gt.construct_gene_tree(True)
# llh_ts = [tree_llh_continous(x) for x in trees_ts]
# print(sum([llh.tree.weight * llh.C_mat() for llh in llh_ts]))

