<!-- README.md is generated from README.Rmd. Please edit that file -->

# phenomix

R package for fitting distributions to run timing data via maximum
likelihood

[![R build
status](https://github.com/ericward-noaa/phenomix/workflows/R-CMD-check/badge.svg)](https://github.com/ericward-noaa/phenomix/actions)

## Installation

You can install phenomix with:

``` r
remotes::install_github("ericward-noaa/phenomix",build_vignettes = TRUE)
```

Load libraries

``` r
library(phenomix)
library(ggplot2)
```

## Functions

The package pheomix provides a suite of curve fitting to describe data
that may be generated from a process when distributions in time might be
concentrated (from fisheries, this occurs with counts over time of
salmon returning from the ocean to spawn or juvenile fish emigrating
from streams to the ocean).

![Predicted (black line) and observed counts (red dots) for hypothetical
dataset. Multiple observations may exist for some days, or no
observations on others.](README-figs/unnamed-chunk-5-1.png)

In a given year, the curve might be described by a symmetric or
asymmetric Gaussian or Student-t distribution (shown here in log-scale
on the y-axis). Questions of interest might be - are the means (x-axis)
shifting through time? - are the variances shifting through time? - does
the model support a symmetric or asymmetric distribution?

![](README-figs/unnamed-chunk-6-1.png)<!-- -->

## Examples

The main functions are `create_data()` and `fit()`. See `?create_data`
and `?fit` for additional details and examples. A vignette includes
additional detail, and examples of several models as well as function
arguments available
[here](https://github.com/ericward-noaa/phenomix/tree/master/vignettes).

## NOAA Disclaimer

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) |
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
