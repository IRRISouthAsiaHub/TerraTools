## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DFBA)

## -----------------------------------------------------------------------------
group1 <- c(12.90, 10.84, 22.67, 10.64, 10.67, 10.79, 13.55, 10.95, 12.19, 12.76, 10.89, 11.02, 14.27, 13.98, 11.52, 13.49, 11.22, 15.07, 15.74, 19.00)

group2 <- c(4.63, 58.64, 5.07, 4.66, 4.13, 3.92, 3.39, 3.57, 3.56, 3.39)

dfba_median_test(E = group1, 
                 C = group2)


## -----------------------------------------------------------------------------
median_test_results<-dfba_median_test(E = group1, 
                                      C = group2)

dfba_binomial(median_test_results$nEabove, 
              median_test_results$nCabove)

## ----fig.width = 7------------------------------------------------------------
plot(dfba_binomial(median_test_results$nEabove, 
                   median_test_results$nCabove))

## -----------------------------------------------------------------------------
dfba_beta_bayes_factor(a_post = median_test_results$a_post, 
                       b_post = median_test_results$b_post, 
                       method = "point", 
                       H0 = 14/15)

