#!/bin/bash

# this script simulate countious trait of a six-taxon tree in seastaR

sptree=~/spinney/tree/6taxa_low.sptree.nwk;
outfile=~/spinney/trait/6taxa_low.csv;
sigma2=1;


module load r
Rscript run_seastaR.R $sptree $outfile $sigma2
