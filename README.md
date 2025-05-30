# Branching rate heterogeneity in lineage tracing trees
This repository contains the essential scripts that support the manuscript titled
'Detecting branching rate heterogeneity with tree balance statistics in lineage tracing trees'.

Authors: Yingnan Gao, Alison F. Feder

DOI: https://doi.org/10.1101/2024.06.27.601073

The accompanying data can be found in repository: TBD.

A breakdown of the scripts can be found below:

## Scripts/Python/
This subfolder holds the .py scripts for generating lineage tracing trees via
Cassiopeia (https://github.com/YosefLab/Cassiopeia). See reproduce_simulation.sh
for more details.

### Scripts/Python/infer_LT_tree_from_CM.py
This script uses Cassiopeia to reconstruct the lineage tracing tree from a character
matrix in a text file using the MaxCut algorithm.

### Scripts/Python/simulate_genealogical_trees.py
This script uses Cassiopeia to simulate genealogical trees of certain sizes and
a growth model, and outputs them as text Newick files.

### Scripts/Python/simulate_lineage_tracing_trees.py
This script uses Cassiopeia to simulate lineage tracing trees from genealogical
trees by simulating character matrices and reconstruct trees from them using the
MaxCut algorithm.

## Scripts/Rscripts/
This subfolder holds the .R scripts for tree generation and other functions to
be run in batch (e.g., on a HPC platform). See Scripts/reproduce_simulation.sh
for more details. The two scripts not used for simulating trees are:

### Scripts/Rscripts/calculate_editing_proportion.R
This script extract the barcode saturation status from the simulated character
matrices. The aggregated output from this script is stored in Simulated_data/stats

###Scripts/Rscripts/infer_for_each_CM_file.R
This script generates the empirical trees using the MaxCut algorithm in the folders
in Empirical_data/ from the character matrix files.

## Scripts/R/
This subfolder holds the .R files to be sourced by other R scripts or in
RStudio. See the comments in each file for more details. 

### Scripts/R/Jq_tree_balance.R
This file contains the functions to calculate J^1 and the Sackin index.

### Scripts/R/phylo_ops.R
This file contains utility functions such as the wrapper function to apply a 
function across a phylogeny.

### Scripts/R/test_of_constant_rate.R
This file contains the function to apply the statistical test against EBR 
using J^1 or the Sackin index.

### Scripts/Analysis_*.R
These files are scripts that reproduces the analyses in the manuscript. They can
be sourced from RStudio with the working directory set to the same directory as
this README.txt.

### Scripts/Fig_*.R
These files are scripts that reproduce the figures in the manuscript.
See Scripts/reproduce_figures_and_tables.sh for more details.

### Scripts/reproduce_figures_and_tables.sh
This is a tutorial to reproduce the figures and tables in the manuscript.

### Scripts/reproduce_simulation.sh
This is a tutorial to reproduce the trees in the manuscript.

