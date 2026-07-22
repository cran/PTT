# PTT: Pólya tree type models

`PTT` provides Bayesian nonparametric density and conditional-density models
built from recursive partitions:

- adaptive Pólya trees (APT) and Markov APTs;
- optional Pólya trees (OPT); and
- conditional APT and conditional OPT models.

The implementation supports univariate and multivariate sample spaces, exact
forward-backward inference on a finite partition tree, posterior predictive
densities, hierarchical MAP partitions, and posterior partition samples.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("MaStatLab/PTT")
```

The package compiles C++ code and therefore requires the standard R build tools
for your platform.

## Marginal density estimation

```r
library(PTT)

set.seed(12345)
x <- c(rbeta(150, 3, 10), rbeta(100, 12, 4))
x_grid <- seq(0.01, 0.99, length.out = 200)

fit <- apt(
  X = x,
  Xpred = x_grid,
  max.resol = 7,
  n.s = 5,
  n.grid = 5
)

plot(x_grid, fit$predictive_densities, type = "l",
     xlab = "x", ylab = "Posterior predictive density")
fit$part_points_hmap
```

Use `opt()` for an optional Pólya tree with the same interface:

```r
opt_fit <- opt(x, Xpred = x_grid, max.resol = 7)
```

## Conditional density estimation

Prediction locations are paired: row `i` of `Xpred` is evaluated with row `i`
of `Ypred`.

```r
set.seed(12345)
x <- runif(300)
y <- ifelse(x < 0.5, rbeta(300, 3, 9), rbeta(300, 9, 3))

y_grid <- seq(0.02, 0.98, length.out = 100)
x_grid <- rep(c(0.25, 0.75), each = length(y_grid))
y_pred <- rep(y_grid, 2)

conditional_fit <- cond.opt(
  X = x,
  Y = y,
  Xpred = x_grid,
  Ypred = y_pred,
  max.resX = 5,
  max.resY = 7
)

matplot(
  y_grid,
  matrix(conditional_fit$predictive_densities, ncol = 2),
  type = "l", lty = 1,
  xlab = "y", ylab = "Conditional predictive density"
)
legend("top", c("x = 0.25", "x = 0.75"), lty = 1, col = 1:2)
```

## Sample spaces and outputs

The default `Omega.type = "unit"` (or `OmegaX.type`/`OmegaY.type`) requires
all corresponding values to lie in `[0, 1]`. Use `"standardized"` to construct
a padded hyperrectangle from the observed and prediction ranges. The fitted
sample-space bounds are returned as `Omega`, or as `OmegaX` and `OmegaY`.

The main returned quantities are:

- `logphi`: log marginal likelihood;
- `logrho`: log posterior complete-shrinkage probability at the root;
- `predictive_densities`: densities in the physical coordinates of the chosen
  sample space;
- `part_points_hmap`: hierarchical MAP partition encoded on the integer grid at
  the requested maximum resolution; and
- `part_points_post_samples`: posterior partition samples, when requested.

Seven installed demos emphasize visual interpretation and use reproducible
simulations. List them with `demo(package = "PTT")`. Three concise gallery
demos are:

```r
demo("marginal_1d_comparison", package = "PTT")
demo("conditional_1d_surface", package = "PTT")
demo("marginal_2d_partition", package = "PTT")
```

The demos show an APT/OPT comparison with terminal-cell resolution, a
conditional-density heatmap with fitted slices, and a bivariate
predictive-density map alongside its hMAP terminal partition.

Four additional mixture and conditional-regime workflows are available with
the same visual style and consistent naming convention:

```r
demo("marginal_1d_mixture", package = "PTT")
demo("marginal_2d_mixture", package = "PTT")
demo("conditional_1d_regimes", package = "PTT")
demo("conditional_2d_mixture", package = "PTT")
```

These retain the original models, simulation sizes, and principal fitting
settings while adding fitted-versus-true comparisons, coordinated palettes,
conditional-density slices, and clearer partition displays.

## References

Ma, L. (2017). Adaptive shrinkage in Pólya tree type models. *Bayesian
Analysis*, 12(3), 779-805. [doi:10.1214/16-BA1021](https://doi.org/10.1214/16-BA1021)

Ma, L. (2017). Recursive partitioning and multi-scale modeling on conditional
densities. *Electronic Journal of Statistics*, 11(1), 1297-1325.
[doi:10.1214/17-EJS1254](https://doi.org/10.1214/17-EJS1254)

Wong, W. H. and Ma, L. (2010). Optional Pólya tree and Bayesian inference.
*The Annals of Statistics*, 38(3), 1433-1459.
[doi:10.1214/09-AOS755](https://doi.org/10.1214/09-AOS755)
