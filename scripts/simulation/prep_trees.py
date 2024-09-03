from ete3 import Tree
import random
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--species_tree', help='file contains species tree (newick format) in coalescent unit', type=str)
    parser.add_argument('--ms_output', help='file contains sampled gene trees genereated by ms', type=str)
    parser.add_argument('--num_gene_trees', help='number of sampled gene trees for inference', type=int)
    parser.add_argument('--output_prefix', help='prefix of output file for sampled gene trees for inference', type=str)

    args = parser.parse_args()
    # Load trees
    sptree_path= args.species_tree
    sptree = Tree(sptree_path, format=1)

    with open(args.ms_output,'r') as f:
        sampled_nwk = f.readlines()
    sampled_nwk = [x for x in sampled_nwk if x.startswith("(")]
    sampled_trees = [ Tree(x, format=1) for x in sampled_nwk]

    # calculate observed discordance
    pairwise_rf = [ sptree.compare(x, unrooted=False)["rf"]  for x in sampled_trees]
    discordance = (1 - sum([True for x in pairwise_rf if x == 0.0])/len(sampled_nwk))*100
    print(f'Observed discordance in ms output : {discordance:.1f}% ({len(sampled_nwk)} gene trees)')

    # sampled gene trees for inference
    output_nwk = random.sample(sampled_nwk, args.num_gene_trees)
    nwk_sub = [ '[&R] %s'%x for x in output_nwk]
    nwk_string = ''.join(nwk_sub)
    # each gene tree has the same weight
    freq_string = '\n'.join(  [str(1/len(nwk_sub)),]*len(nwk_sub)  )

    pairwise_rf = [ sptree.compare(Tree(x), unrooted=False)["rf"]  for x in output_nwk]
    discordance = (1 - sum([True for x in pairwise_rf if x == 0.0])/len(output_nwk))*100
    print(f'Observed discordance in {len(nwk_sub)} sampled gene trees: {discordance:.1f}%')

    with open('%s.sampled.txt'%args.output_prefix, 'w') as f:
        template = ['begin genetrees;', nwk_string, 'end;',
                    'begin genetreefreqs;', freq_string, 'end;']
        f.write('\n'.join(template))


    with open('%s.sptree.txt'%args.output_prefix, 'w') as f:
        template = ['begin genetrees;', '[&R] ' + sptree.write(), 'end;',
                    'begin genetreefreqs;', str(1.0), 'end;']
        f.write('\n'.join(template))


    with open('%s.ratetree_1p.nwk'%args.output_prefix, 'w') as f:
        for node in sptree.traverse():
            node.dist = 1
        template = ['[&R]', sptree.write()]
        f.write(' '.join(template))