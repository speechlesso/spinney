#!/bin/bash

#SBATCH -J expectandsp_1p
#SBATCH -p general
#SBATCH -o high_2_1p.txt
#SBATCH -e high_2_1p.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=moyu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH -A r00526

export OMP_NUM_THREADS=4

condition=high_2;

traitfile=6taxa_${condition}.csv;
gtree=6taxa_${condition}.expected.txt;
sptree=6taxa_${condition}.sptree.nwk;
cdtree=6taxa_${condition}.concord.nwk;
numtree=4;

cd ~/spinney
module load python/3.9

### gene trees
#1-param
python ./code/sampledtree_rates.py -f ./trait/$traitfile -N $numtree -genetrees ./tree/$gtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_1p.ratetree.nwk -condition ${condition}_expected_1p

#2-param
#python ./code/sampledtree_rates.py -f ./trait/$traitfile -N $numtree -genetrees ./tree/$gtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_2p.ratetree.nwk -condition ${condition}_expected_2p

### species tree
python ./code/sampledtree_rates.py -f ./trait/$traitfile -N 1 -genetrees ./tree/$sptree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_1p.ratetree.nwk  -condition ${condition}_sp_1p

#python ./code/sampledtree_rates.py -f ./trait/$traitfile -N 1 -genetrees ./tree/$sptree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_2p.ratetree.nwk  -condition ${condition}_sp_2p

### concordant tree
#python ./code/sampledtree_rates.py -f ./trait/$traitfile -N 1 -genetrees ./tree/$cdtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_1p.ratetree.nwk -condition ${condition}_concord_1p

#python ./code/sampledtree_rates.py -f ./trait/$traitfile -N 1 -genetrees ./tree/$cdtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_2p.ratetree.nwk -condition ${condition}_concord_2p
