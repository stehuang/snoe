

# Nonlinear Causal Learning Algorithm

R package for a causal discovery method focused on learning non-linear relations from data

Author: Stella Huang (stellahyh@g.ucla.edu)

## Installation

Install this package from Github by running the following code in R:

```r
devtools::install_github("stehuang/nleo")
```

## Summary

- Highlights: 

	- A likelihood-ratio based statistical test to compare and determine the true causal direction of an undirected edge;
	- A framework for ranking undirected edges for orientation based on an independence measure derived from the Additive Noise Model assumption;
	- Faster computation times and greater interpretability of results during the learning process compared to competing non-linear learning methods;
	- Theoretical guarantees for structural learning in both the population and finite sample settings.


- Algorithm Structure:

	1. Initial structural learning: performs modified version of PC algorithm to learn a denser CPDAG containing extraneous edges but v-structures/edge orientations detected under a strigent threshold
	2. Edge orientation: employs likelihood-ratio based statistical test to determine the true orientation of undirected edges in the graph
	3. Edge deletion: utilizes significance test results of covariates to delete superfluous edges from the graph, yielding the final DAG structure
	4. (For practical purposes) Graph refinement: converts previous output from PDAG to DAG by applying the likelihood test to undirected edges, if there exists any


## Example

## Dependencies

This algorithm utilizes functions from various packages, with the main ones being:

- package `bnlearn`: modifications made by building upon the PC algorithm, finding conditional independence relations, learning edge orientations, and evaluating learned graphs
- package `mgcv`: implemented Generalized Additive Models (GAM) as our choice of the regression method
- package `infotheo`: discretizes the data and calculated the pairwise mutual information between two variables for the conditional independence test in the PC algorithm







