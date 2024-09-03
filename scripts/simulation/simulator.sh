#!/bin/bash

sptreenwk="((((A0:3.7,A1:3.7):0.1,(A2:3.7,A3:3.7):0.1):0.1,A4:3.9):0.1,A5:4);";
sigma2=1;
total_gene_tree=1000;
sampled_gene_tree=20;
output_prefix=temp
outfile=./${output_prefix}.csv;

# Read the template and replace the placeholders with actual values
sed "s/SPTREE/$sptreenwk/g; s/NUM_GENE_TREE/$total_gene_tree/g"\
 run_phylonet_template.nex > run_phylonet.nex

echo "############## Simulation settings"
echo "sigma2:" $sigma2
echo "species tree:" $sptreenwk
echo "total number of sampled gene trees:" $total_gene_tree
echo "number of trees for inference:" $sampled_gene_tree
echo "output trait file:" $outfile
echo $sptreenwk > "sptree.nwk"

# simulate sufficient gene trees with ms
java -jar ~/bin/PhyloNet.jar run_phylonet.nex > ${output_prefix}.ms
# add prefix for each line
#sed -i 's/^/[\&R] /' temp.ms
echo

echo "############## Gene trees info"
python prep_trees.py --species_tree sptree.nwk --ms_output ${output_prefix}.ms \
--num_gene_trees $sampled_gene_tree --output_prefix $output_prefix
echo

# simulate traits with seastaR
echo "############## seastaR inference"
module load r
Rscript run_seastaR.R sptree.nwk $outfile $sigma2

# Replace prefix "abc" with "xyz"
#for file in temp.*; do mv "$file" "highILS_condition7${file#temp}"; done