#!/usr/bin/env python3
"""
This script generates ancestral state maps for a single trait
using the example data included in this repository.

This script generates ancestral states for the 1st trait with
the species tree and Spinney respectively.

Usage:
    python ancestral_state.py
"""


import matplotlib.pyplot as plt
from utils import ancestral_state_pipeline


CONFIG = {"condition": "example",  
          "speciation_time": 1.25,  
          "final_time": 3, 
          "num_slicing": 100, 
          "traits": [0], 
          # Change this to absolute directory to "example" folder
          "base_dir": "" 
          }


fig, axes = plt.subplots(1, 2, figsize=(20, 8), dpi=100)

_, cs1 = ancestral_state_pipeline(
    CONFIG,
    prefix="example",
    treetype="st",
    trait_ids=[0],
    ax=axes[0],
)

_, cs2 = ancestral_state_pipeline(
    CONFIG,
    prefix="example",
    treetype="sy",
    trait_ids=[0],
    ax=axes[1],
)

plt.savefig("example_ancestral_state_map.png")