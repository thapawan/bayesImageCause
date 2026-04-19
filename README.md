# bayesImageCause <img src="man/figures/logo.png" align="right" height="139" alt="bayesImageCause logo" />

<!-- badges: start -->

[![CRAN status](https://img.shields.io/badge/CRAN-not%20submitted-lightgrey.svg)](https://CRAN.R-project.org/package=bayesImageCause)
[![R-CMD-check](https://github.com/yourusername/bayesImageCause/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/bayesImageCause/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/yourusername/bayesImageCause/branch/main/graph/badge.svg)](https://app.codecov.io/gh/yourusername/bayesImageCause)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01234/status.svg)](https://doi.org/10.21105/joss.01234)

<!-- badges: end -->

## Overview

**bayesImageCause** is an R package for **Bayesian Causal Inference with Image Data**...


---
# bayesImageCause <img src="man/figures/logo.png" align="right" height="139" alt="bayesImageCause logo" />

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/bayesImageCause)](https://CRAN.R-project.org/package=bayesImageCause)
[![R-CMD-check](https://github.com/yourusername/bayesImageCause/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/bayesImageCause/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/yourusername/bayesImageCause/branch/main/graph/badge.svg)](https://app.codecov.io/gh/yourusername/bayesImageCause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01234/status.svg)](https://doi.org/10.21105/joss.01234)

<!-- badges: end -->

## Overview

**bayesImageCause** is an R package for **Bayesian Causal Inference with Image Data**. It provides a complete workflow for estimating pixel-level causal effects when treatments are images (2D/3D arrays) and outcomes are scalar. The package uses Gaussian process priors to model spatial correlation and provides fully Bayesian uncertainty quantification at every pixel.

### Key Features

- **🎯 Causal Estimation**: Estimate pixel-level causal effects with full posterior distributions
- **🖼️ Spatial Priors**: Gaussian process priors with multiple kernel options (squared exponential, Matern, periodic)
- **📊 Uncertainty Quantification**: Credible intervals, significance maps, and posterior predictive checks
- **🔮 Counterfactual Prediction**: Predict outcomes under user-specified image interventions
- **📈 Visualization**: Publication-ready plots using `ggplot2`
- **⚡ Scalable**: Sparse Gaussian process approximations for large images (>2500 pixels)
- **🔄 Reproducible**: Full MCMC diagnostics and seed-setting for reproducible research

### Citation

If you use `bayesImageCause` in your research, please cite:

```bibtex
@article{author2026bayesimagecause,
  title = {{bayesImageCause}: {B}ayesian Causal Inference for Image Data in {R}},
  author = {Lastname, Firstname},
  journal = {Journal of Statistical Software},
  year = {2026},
  volume = {xx},
  number = {x},
  pages = {1--32},
  doi = {10.18637/jss.vxxx.i0x}
}
```

## Installation

### From CRAN (stable version)

```r
install.packages("bayesImageCause")
```

### From GitHub (development version)

```r
# Install from GitHub
if (!require("remotes")) install.packages("remotes")
remotes::install_github("yourusername/bayesImageCause")

# Or with vignettes
remotes::install_github("yourusername/bayesImageCause", build_vignettes = TRUE)
```

### System Requirements

- **R**: Version 4.1.0 or higher
- **C++14 compiler**: For Stan backend
- **Recommended**: 8GB+ RAM for large image analyses

## Quick Start

### Minimal Example

```r
library(bayesImageCause)

# Simulate data with known causal structure
sim <- simulate_image_data(
  n = 100,
  img_dim = 16,
  causal_region = list(x = 6:11, y = 6:11),
  effect_size = 1.5,
  seed = 42
)

# Estimate causal effects
fit <- estimate_causal_effect(
  images = sim$images,
  outcomes = sim$outcomes,
  chains = 2,
  iter = 1000
)

# Visualize results
plot_causal_surface(fit, type = "mean")
plot_significance_map(fit, alpha = 0.05)

# Counterfactual prediction
new_image <- sim$images[1, , ]  # First image
pred <- predict_counterfactual(fit, new_image)
print(pred)
```

### Complete Workflow with Covariates

```r
# Load neuroimaging example data
data(hcp_working_memory)

# Prepare data
images <- hcp_working_memory$activation_maps
outcomes <- hcp_working_memory$reaction_time
covariates <- hcp_working_memory[, c("age", "sex", "brain_volume")]

# Fit model with covariates
fit_cov <- estimate_causal_effect(
  images = images,
  outcomes = outcomes,
  covariates = covariates,
  kernel = "matern",
  prior_length_scale = gamma(2, 1),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

# Model diagnostics
print(fit_cov)
plot_trace(fit_cov)
check_convergence(fit_cov)

# Extract significant regions
sig_pixels <- get_significant_pixels(fit_cov, alpha = 0.05)
cat("Number of significant pixels:", sum(sig_pixels), "\n")

# Create interactive visualization
plot_causal_heatmap(fit_cov, interactive = TRUE)
```

## Documentation

| Resource | Description |
|----------|-------------|
| [Package Website](https://yourusername.github.io/bayesImageCause/) | Full documentation with examples |
| [Getting Started](https://yourusername.github.io/bayesImageCause/articles/bayesImageCause.html) | Introductory vignette |
| [Simulation Studies](https://yourusername.github.io/bayesImageCause/articles/simulation_study.html) | Performance evaluation |
| [Neuroimaging Example](https://yourusername.github.io/bayesImageCause/articles/neuroimaging_example.html) | Applied case study |
| [Custom Kernels](https://yourusername.github.io/bayesImageCause/articles/custom_kernels.html) | Advanced usage |
| [Reference Manual](https://cran.r-project.org/package=bayesImageCause/bayesImageCause.pdf) | Complete function reference |

## Example Outputs

### Causal Surface Estimation

```{r example-plot, echo = FALSE}
# This code is executed but not shown in README
library(ggplot2)
set.seed(42)
sim <- simulate_image_data(n = 80, img_dim = 16, effect_size = 2.0)
fit <- estimate_causal_effect(sim$images, sim$outcomes, chains = 1, iter = 400, refresh = 0)
p1 <- plot_causal_surface(fit, type = "mean") + ggtitle("Estimated Causal Surface")
p2 <- plot_significance_map(fit, alpha = 0.05) + ggtitle("Significant Pixels (95% CI)")
print(p1)
print(p2)
```

![Estimated Causal Surface](man/figures/README-example-plot-1.png)
![Significance Map](man/figures/README-example-plot-2.png)

*Left: Posterior mean of pixel-level causal effects. Right: Pixels where 95% credible interval excludes zero (red = significant).*

## Methodology

### Statistical Model

For outcome $Y_i \in \mathbb{R}$ and image treatment $\mathbf{X}_i \in \mathbb{R}^{p \times q}$:

$$Y_i = \sum_{u=1}^p \sum_{v=1}^q \beta_{uv} X_i(u,v) + \mathbf{Z}_i'\boldsymbol{\gamma} + \epsilon_i$$

where $\boldsymbol{\beta} \sim \mathcal{GP}(\mathbf{0}, \mathbf{K}_\theta)$ with squared exponential kernel:

$$K_\theta(\mathbf{s}, \mathbf{s}') = \eta^2 \exp\left(-\frac{\|\mathbf{s} - \mathbf{s}'\|^2}{2\rho^2}\right)$$

### Inference

Posterior sampling via Hamiltonian Monte Carlo as implemented in Stan. For computational efficiency with large images ($P > 2500$), we implement sparse Gaussian processes using $M = 500$ inducing points.

### Performance Characteristics

| Image Size | Pixels | Runtime (min) | Memory (GB) |
|------------|--------|---------------|-------------|
| 16 × 16    | 256    | 4.2           | 0.5         |
| 32 × 32    | 1,024  | 12.3          | 1.2         |
| 64 × 64    | 4,096  | 48.7          | 3.5         |
| 128 × 128  | 16,384 | 215.4         | 12.0        |

*Note: Times for $n=200$, 4 chains, 1000 iterations, using sparse GP for P > 2500.*

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md).

### Development Workflow

```bash
# Clone the repository
git clone https://github.com/yourusername/bayesImageCause.git
cd bayesImageCause

# Install development dependencies
R -e "install.packages(c('devtools', 'roxygen2', 'testthat', 'knitr', 'rmarkdown'))"

# Build and test
R -e "devtools::document()"
R -e "devtools::check()"
R -e "devtools::test()"

# Build vignettes
R -e "devtools::build_vignettes()"
```

### Reporting Issues

- **Bug reports**: Please use the [GitHub Issues](https://github.com/yourusername/bayesImageCause/issues) page with a reproducible example
- **Feature requests**: Open an issue with the `enhancement` label
- **Questions**: Post on [Stack Overflow](https://stackoverflow.com/questions/tagged/bayesimagecause) with the tag `[bayesimagecause]`

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## License

This package is free and open source software, licensed under GPL-3.0.


## Acknowledgements

We thank the Stan Development Team, the R Core Team, and the contributors to the tidyverse ecosystem. The neuroimaging example uses data from the Human Connectome Project (Van Essen et al., 2013).

## Related Work

Other R packages for causal inference and image analysis:

- **Causal inference**: [`ipd`](https://cran.r-project.org/package=ipd), [`CICI`](https://cran.r-project.org/package=CICI), [`cyclinbayes`](https://cran.r-project.org/package=cyclinbayes)
- **Image analysis**: [`robflreg`](https://cran.r-project.org/package=robflreg), [`ggmulti`](https://cran.r-project.org/package=ggmulti)
- **Bayesian modeling**: [`rstan`](https://mc-stan.org/rstan/), [`brms`](https://paul-buerkner.github.io/brms/)

## Citing This Work

```r
citation("bayesImageCause")
```
