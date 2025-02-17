

# Sequential Nonlinear Edge Orientation (SNEO)

R package for a causal discovery method focused on learning non-linear relations from data

Author: Stella Huang (stellahyh@g.ucla.edu)

## Installation

Install this package from Github by running the following code in R:

```r
devtools::install_github("stehuang/sneo")
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

Below is an example to generate data given a DAG, learn a DAG using the algorithm, and evaluating the results

```{r}
# generate data
true.dag <- get_true_dag("asia")
data <- DAGdatagen(n = 1000, dag_amat=amat(true.dag), data_fx='inv', ninv.method="random", se=0.5)$DAGdata

# run algorithm
# `ll_test_approach` options are sample-splitting (SS) or cross-validation (CV)
adjmat <- sneo(data=data, ll_test_approach="SS")

# compare with ground truth
learning_res <- get_f1(amat(true.dag), list(adjmat))
names(learning_res) <- c("shd", "f1", "tp", "fp", "fn", "wrong_dir", "wrong_dir_undir", "n_undir")
print(learning_res)
```

## Dependencies

This algorithm utilizes functions from various packages, with the main ones being:

- package `bnlearn`: modifications made by building upon the PC algorithm, finding conditional independence relations, learning edge orientations, and evaluating learned graphs
- package `mgcv`: implemented Generalized Additive Models (GAM) as our choice of the regression method
- package `infotheo`: discretizes the data and calculated the pairwise mutual information between two variables for the conditional independence test in the PC algorithm







