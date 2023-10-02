#!/bin/bash

#SBATCH -J 6taxa_low
#SBATCH -p general
#SBATCH -o 6taxa_low.txt
#SBATCH -e 6taxa_low.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=1G
#SBATCH -A general

condition=low;

traitfile=6taxa_${condition}.csv;
gtree=6taxa_${condition}.ms;
sptree=6taxa_${condition}.sptree.nwk;
numtree=20;

cd ~/spinney
#module load python/3.9

#1-param
#python3 ./code/sampledtree_rates.py -f ./trait/$traitfile -N $numtree -genetrees ./tree/$gtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_1p.ratetree.nwk -condition {$condition}_1p

#2-param
python3 ./code/sampledtree_rates.py -f ./trait/$traitfile -N $numtree -genetrees ./tree/$gtree -speciestree ./tree/$sptree -ratetree ./tree/6taxa_2p.ratetree.nwk -condition ${condition}_2p
