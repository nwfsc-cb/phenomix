## ----setup, include = FALSE, cache=FALSE--------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.618
)

## ----packages, message=FALSE, warning=TRUE------------------------------------
library(ggplot2)
library(salmix)
library(dplyr)
library(TMB)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
set.seed(123)
df <- data.frame(x = seq(110, 150, 1))
df$y <- dnorm(df$x, mean = 130, 10) * 10000
df$obs <- rnorm(nrow(df), df$y, 30)
ggplot(df, aes(x, obs)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  xlab("Day of year") +
  ylab("Observation")

## -----------------------------------------------------------------------------
glimpse(fishdist)

## -----------------------------------------------------------------------------
names(fishdist)

## -----------------------------------------------------------------------------
datalist <- create_data(fishdist,
  min_number = 0,
  variable = "number",
  time = "year",
  date = "doy",
  asymmetric_model = FALSE,
  est_sigma_trend = TRUE,
  est_mu_trend = TRUE,
  est_sigma_re = TRUE,
  est_mu_re = TRUE,
  tail_model = "gaussian",
  family = "gaussian"
)

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
set.seed(1)
fitted <- fit(datalist)

## -----------------------------------------------------------------------------
names(fitted)

## -----------------------------------------------------------------------------
fitted$pars$convergence

## -----------------------------------------------------------------------------
sdrep_df <- data.frame(
  "par" = names(fitted$sdreport$value),
  "value" = fitted$sdreport$value, "sd" = fitted$sdreport$sd
)
head(sdrep_df)

## ----eval = FALSE-------------------------------------------------------------
#  TMB::sdreport(fitted$obj)

## ---- fig.cap="Fitted symmetric model with tails from a  Gaussian distribution", fig.width = 8----
g <- plot_diagnostics(fitted, type = "timing", logspace = TRUE)
g

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
set.seed(2)

datalist <- create_data(fishdist,
  min_number = 0,
  variable = "number",
  time = "year",
  date = "doy",
  asymmetric_model = TRUE,
  tail_model = "student_t"
)
fitted_t <- fit(datalist)

## ---- fig.cap="Fitted asymmetric model with heavy tails from a t-distribution", fig.width = 8----
plot_diagnostics(fitted_t)

## -----------------------------------------------------------------------------
aic_1 <- extractAIC(fitted)$AIC
aic_1

aic_2 <- extractAIC(fitted_t)$AIC
aic_2

## ----message=FALSE, warning=FALSE, results='hide', eval=FALSE-----------------
#  set.seed(5)
#
#  datalist = create_data(fishdist,
#    min_number=0,
#    variable = "number",
#    time="year",
#    date = "doy",
#    asymmetric_model = FALSE,
#    tail_model = "gnorm")
#  fitted = fit(datalist, limits = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  fit = fitted(..., control=list(rel.tol = 1.0e-12,
#    eval.max=4000, iter.max=4000))
