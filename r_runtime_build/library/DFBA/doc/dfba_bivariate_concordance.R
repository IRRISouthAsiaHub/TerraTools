## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DFBA)

## -----------------------------------------------------------------------------
x <- c(47, 39, 47, 42, 44, 46, 39, 37, 29, 42, 54, 33, 44, 31, 28, 49, 32, 37, 46, 55, 31)

y <- c(36, 40, 49, 45, 30, 38, 39, 44, 27, 48, 49, 51, 27, 36, 30, 44, 42, 41, 35, 49, 33)

A <- dfba_bivariate_concordance(x,
                                y)

A

## ----fig.width = 7, fig.height = 4--------------------------------------------
plot(A)

## -----------------------------------------------------------------------------
# predicted values from model

p = seq(.05, .95, .05)

ypred= 17.332 - (50.261*p) + (48.308*p^2)

# Note: the coefficients in the ypred equation were found first via a polynomial regression

# observed values

yobs <- c(19.805, 10.105, 9.396, 8.219, 6.110, 4.543, 5.864, 4.861, 6.136, 5.789,
          5.443, 5.548, 4.746, 6.484, 6.185, 6.202, 9.804, 9.332, 14.408)

B <- dfba_bivariate_concordance(x = yobs, 
                                y = ypred, 
                                fitting.parameters = 3)

B

## ----fig.width = 7, fig.height = 4--------------------------------------------
plot(B)

