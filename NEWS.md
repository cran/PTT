# PTT 1.0.1

## Correctness

- Fixes undefined behavior in conditional hierarchical MAP partition output
  caused by using an unsigned counter with `-1` as its initial sentinel.
- Uses 64-bit node offsets to avoid truncation when indexing larger trees.
- Rejects predictor/response dimension and shrinkage-state configurations that
  exceed the capacity of the internal partition representation.

## Validation

- Adds regression checks for conditional partition boundaries and supported
  input capacities.
- Audits the remaining bit shifts and verifies the package examples, tests,
  and demos with undefined-behavior and address sanitizers.

# PTT 1.0

PTT 1.0 is the first official release of the package.

## Models

- Provides adaptive Pólya tree (APT), Markov APT, and optional Pólya tree
  (OPT) models for Bayesian nonparametric density estimation.
- Provides conditional APT and conditional OPT models for conditional-density
  estimation through recursive partitioning of the predictor space.
- Supports univariate and multivariate observations on the unit hypercube or a
  data-adaptive standardized hyperrectangle.

## Inference and results

- Uses exact forward-backward inference on finite recursive partition trees.
- Computes log marginal likelihoods, posterior root shrinkage probabilities,
  posterior predictive densities, and hierarchical MAP partitions.
- Generates posterior partition samples using R's random-number generator for
  reproducible simulation.
- Reports predictive densities in the physical coordinates of the selected
  sample space.

## Package interface

- Supplies documented R interfaces backed by registered C++ routines using
  Rcpp and RcppArmadillo.
- Validates observations, prediction matrices, sample spaces, and model
  parameters before native computation.
- Includes plotting utilities for two-dimensional partitions and ROC curves.
- Includes package-level documentation, method references, citation metadata,
  API and numerical regression tests, and seven reproducible visual demos for
  marginal and conditional models in one and two dimensions.
