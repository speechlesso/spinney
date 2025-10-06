#!/bin/bash

BASE_DIR="$(git rev-parse --show-toplevel)"
INFERENCE_SCRIPT=${BASE_DIR}/src/spinney/rate_inference.py

cd $BASE_DIR/example

condition="example"
traitfile=${condition}.csv;
genetrees=${condition}.sy.txt;
speciestree=${condition}.st.txt;
ratetree_1p=${condition}.ratetree_1p.nwk
ratetree_2p=${condition}.ratetree_2p.nwk

### species tree
numprocess=1;
python $INFERENCE_SCRIPT -f $traitfile --genetrees $speciestree --speciestree $speciestree --ratetree $ratetree_1p --condition ${condition}_st_1p --num_process $numprocess

### Spinney
numprocess=5;
python $INFERENCE_SCRIPT -f $traitfile --genetrees $genetrees --speciestree $speciestree --ratetree $ratetree_1p --condition ${condition}_sy_1p --num_process $numprocess
