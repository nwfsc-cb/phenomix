---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# salmix
R package for fitting distributions to run timing data via maximum likelihood

[![Build Status](https://travis-ci.com/ericward-noaa/salmix.svg?branch=master)](https://travis-ci.com/ericward-noaa/salmix)

## Installation

You can install salmix with:

``` r
devtools::install_github("ericward-noaa/salmix")
```
## Functions

The package salmix provides a suite of curve fitting to describe data that may be generated from a process when distributions in time might be concentrated (from fisheries, this occurs with counts over time of salmon returning from the ocean to spawn or juvenile fish emigrating from streams to the ocean). 

In a given year, the curve might be described by a symmetric or asymmetric Gaussian or Student-t distribution (shown here in log-scale on the y-axis). Questions of interest might be
- are the means (x-axis) shifting through time?
- are the variances shifting through time?
- does the model support a symmetric or asymmetric distribution?

```{r echo=FALSE}
set.seed(123)
df = expand.grid(year = 1:10, doy=100:250)
mus = rnorm(10, 175, 20)
thetas = runif(10, 5, 8)
sigmas = rnorm(10, 10, 5)
df$y = dnorm(df$doy, mus[df$year], sigmas[df$year], log=TRUE) + thetas[df$year]
df$year = as.factor(df$year)
```

```{r}
ggplot2::ggplot(df, aes(doy,y,col=year,group=year)) + geom_line()
```

## Examples

The main functions are `create_data()` and `fit()`. See `?create_data` and `?fit` for additional details and examples. A vignette includes additional detail, and examples of several models as well as function arguments available [here](https://github.com/ericward-noaa/salmix/tree/master/vignettes).
