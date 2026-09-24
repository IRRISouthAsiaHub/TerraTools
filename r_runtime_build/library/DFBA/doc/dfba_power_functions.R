## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ----setup--------------------------------------------------------------------
library(DFBA)

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(1)
#  A <- dfba_bayes_vs_t_power(delta = 0.3,
#                             model = "normal",
#                             design = "independent",
#                             hide_progress = TRUE)
#  
#  A

## ----echo = FALSE-------------------------------------------------------------
load("power_ex1")
A

## ----fig.width = 7------------------------------------------------------------
plot(A)

## ----eval = FALSE-------------------------------------------------------------
#  B <- dfba_power_curve(n = 40,
#                        model = "normal",
#                        design = "paired",
#                        hide_progress = TRUE)
#  
#  B

## ----echo = FALSE-------------------------------------------------------------
load("power_ex2")
B

## ----fig.width = 7------------------------------------------------------------
plot(B)

## ----eval = FALSE-------------------------------------------------------------
#  C<- dfba_power_curve(n = 40,
#                       model = "weibull",
#                       design = "paired",
#                       shape1=.8,
#                       shape2=.8,
#                       hide_progress = TRUE)
#  
#  C

## ----echo = FALSE-------------------------------------------------------------
load("power_ex3")
C

## ----fig.width=7--------------------------------------------------------------
plot(C)

