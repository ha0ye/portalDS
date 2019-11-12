Rocky Intertidal Dynamic Stability Analysis
================
Hao Ye
2019-11-12

# Introduction

This report documents applying the dynamic stability analysis to the
rocky intertidal system reported on in (Benincà et al. 2015)

First, some setup for the environment:

``` r
library(portalDS)
library(dplyr)
library(ggplot2)

set.seed(42)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Workflow

## Data

We use data from the 8 control plots, up to the 2015 treatment-switch.
Abundances are summed, and scaled according to the number of plots
sampled in each census (sometimes, incomplete sampling occurred). Only
species with at least 50% non-zero abundance across the time series were
kept. Missing censuses had abundances imputed using linear interpolation
for each individual species.

``` r
block <- make_portal_block(filter_q = 0.5)
str(block)
```

## Analysis

We do not go through the full analysis here. Instead, see the [Maizuru
Dynamic Stability
vignette](https://ha0ye.github.io/portalDS/articles/maizuru-dynamic-stability.html)
or (eventual methdos write-up).

We specify a results file for storing the outputs of the analysis. The
`compute_dynamic_stability()` function will check for the presence of
individual output components and will skip them if they’ve already been
computed.

``` r
results_file <- here::here("output/portal_ds_results_50.RDS")
results <- compute_dynamic_stability(block, results_file)
str(results, max.level = 1)
```

# Results

## Abundance time series

``` r
plot_time_series(results$block)
```

## Eigenvalues & Eigenvectors

Highlighted segments are the posterior estimates of regime shifts from
Christensen et al. 2018.

``` r
plot_eigenvalues(results$eigenvalues, 
                 num_values = 3) %>% 
    add_regime_shift_highlight()
```

``` r
plot_eigenvectors(results$eigenvectors) %>% 
    add_regime_shift_highlight()
```

## Singular values and SVD vectors

Again, highlighted segments are the posterior estimates of regime shifts
from Christensen et al. 2018.

``` r
plot_svd_values(results$svd_decomp$d, 
                num_values = 3) %>% 
    add_regime_shift_highlight()
```

``` r
plot_svd_vectors(results$svd_decomp$u) %>% 
    add_regime_shift_highlight()
```

## Volume contraction and total variance

``` r
plot_volume_contraction(results$volume_contraction) %>%
  add_regime_shift_highlight()
```

``` r
plot_total_variance(results$total_variance) %>%
  add_regime_shift_highlight()
```

<div id="refs" class="references">

<div id="ref-Beninca_2015">

Benincà, Elisa, Bill Ballantine, Stephen P. Ellner, and Jef Huisman.
2015. “Species Fluctuations Sustained by a Cyclic Succession at the Edge
of Chaos.” *Proceedings of the National Academy of Sciences* 112 (20).
Proceedings of the National Academy of Sciences: 6389–94.
<https://doi.org/10.1073/pnas.1421968112>.

</div>

</div>
