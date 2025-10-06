# Spinney
Spinney is a Python package for inferring evolutionary rate variation across species while accounting for gene tree discordance. It implements a Brownian motion–based model that integrates information from multiple gene trees, species trees, and observed continuous traits. Spinney can estimate branch-specific rates and reconstruct ancestral states, providing a robust framework for studying trait evolution under incomplete lineage sorting (ILS).

# Install dependency
```
git clone https://github.com/speechlesso/spinney.git
cd spinney
pip install -r requirements.txt
```
# Rate inference
## Input

### Trait file
A CSV file where the header contains species names, and each row corresponds to one trait.
### Trees
- Species tree  
A text file containing a rooted species tree in newick, whose weight is 1.  

- Gene trees  
A text file containing multiple rooted gene trees, each with an associated weight.  

- Rate tree   
A newick tree with the same topology as the species tree. The branch lengths represent the rate assignments for inference.  

## How to Run
Run module `./src/spinney/rate_inference.py`. An example [shell script](scripts/inference.sh) demonstrates how to specify input files and parameters.

## Output
The output is saved as a JSON file summarizing the optimization results.
```
{
  "n_trees": 20,
  "minimize_success": true,
  "rate_0": 3.3593460365720733,
  "condition": "example_sy_1p",
  "trait_idx": 0
}
```

# Ancestral state inference
After obtaining rate estimates, we can reconstruct ancestral state at given speciation event [(python script)](src/spiNney/ancestral_state.py).




