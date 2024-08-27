#!/bin/bash

#SBATCH -J 6taxa_cond1_gt
#SBATCH -p general
#SBATCH -o 6taxa_cond1_gt_%j.txt
#SBATCH -e 6taxa_cond1_gt_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=moyu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=00:20:00
#SBATCH --mem=1G
#SBATCH -A r00526

#Load any modules that your program needs
source ~/spinney/spinney_env/bin/activate

#same as --ntasks-per-node
numprocess=5;

for condition in "6taxa_condition1";
do
    cd ~/spinney/experiment/$condition
    traitfile=${condition}.csv;
    gtree=${condition}.sampled.txt;
    sptree=${condition}.sptree.txt;
    cdtree=${condition}.concord.txt;
    numtree=10;

    echo "" > ${condition}.tsv;
    pypath=../../code/sampledtree_rates.py

    ### gene trees
    #1-param
    python $pypath -f ./$traitfile -N $numtree -genetrees ./$gtree -speciestree ./$sptree -ratetree ../6taxa_1p.ratetree.nwk \
    -condition ${condition}_sampled_1p -num_process $numprocess

    #2-param
    python $pypath -f ./$traitfile -N $numtree -genetrees ./$gtree -speciestree ./$sptree -ratetree ../6taxa_2p.ratetree.nwk \
    -condition ${condition}_sampled_2p -num_process $numprocess

done
