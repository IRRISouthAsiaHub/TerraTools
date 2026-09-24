library(dplyr)
library(survival)
set.seed(1)

make_weights <- function(x){
    z <- runif(nrow(x))
    z / mean(z)
}

pprint_check <- function(mod1, mod2, name = NULL){
    res <- vector(mode = "list", length = 7)
    names(res) <- c("coefficients", "deviance_residuals", "pearson_residuals",
                    "response_residuals", "fitted_values", "deviance", "log_likelihood")
    diffs <- function(one, two) c(abs = max(abs(one - two)), rel = max(abs( (one - two) / (one + two) / 2)))
    res$coefficients <- diffs(coef(mod1, subset = "covariates"), coef(mod2))
    res$deviance_residuals <- diffs(resid(mod1, type = "deviance"), resid(mod2, type = "deviance"))
    res$pearson_residuals <- diffs(resid(mod1, type = "pearson"), resid(mod2, type = "pearson"))
    res$response_residuals <- diffs(resid(mod1, type = "response"), resid(mod2, type = "response"))
    res$fitted_values <- diffs(fitted(mod1), fitted(mod2))
    res$deviance <- diffs(deviance(mod1), deviance(mod2))
    res$log_likelihood <- diffs(logLik(mod1), logLik(mod2))
    if (! is.null(name)) cat("\n- ", name, "\n\n", sep = "")
    stars <- function(x) case_when(x < 1E-06             ~ ".",
                                   x > 1E-06 & x < 1E-05 ~ "*",
                                   x > 1E-05 & x < 1E-04 ~ "**",
                                   x > 1E-04 & x < 1E-03 ~ "***",
                                   x > 1E-03 & x < 1E-02 ~ "****",
                                   x > 1E-02 & x < 1E-01 ~ "*****",
                                   x > 1E-01              ~ "******")
    for (i in 1:length(res)){
        cat(format(names(res)[i], width = 20), ":",
            "abs: " , format(res[[i]]["abs"], digits = 1, scientific = TRUE),
            " rel: ", format(res[[i]]["rel"], digits = 1, scientific = TRUE),
            " ", stars(res[[i]]["rel"]),
                        "\n", sep = "")
    }
    invisible(res)
}    

check_glm <- function(mod1, mod2, tolerance = 1E-04, name, skip = NULL){
    .tolerance <- tolerance
    .skip <- rep(TRUE, 8)
    if (! is.null(skip)) .skip[skip] <- FALSE
    compar <- intersect(names(coef(mod1)), names(coef(mod2)))
    if(.skip[1]) expect_equal(coef(mod1)[compar], coef(mod2)[compar], tolerance = .tolerance,
                              info = paste(name, "coefficients", sep = ": "))
    if(.skip[2]) expect_equal(resid(mod1), resid(mod2), tolerance = .tolerance,
                               info = paste(name, "residuals", sep = ": "))
    if(.skip[3]) expect_equal(resid(mod1, type = "deviance"), resid(mod2, type = "deviance"), tolerance = .tolerance,
                              info = paste(name, "deviance residuals", sep = ": "))
    if(.skip[4]) expect_equal(resid(mod1, type = "pearson"), resid(mod2, type = "pearson"), tolerance = .tolerance,
                              info = paste(name, "pearson residuals", sep = ": "))
    if(.skip[5]) expect_equal(resid(mod1, type = "response"),  resid(mod2, type = "response"), tolerance = .tolerance,
                              info = paste(name, "response residuals", sep = ": "))
    if(.skip[6]) expect_equal(fitted(mod1), fitted(mod2), tolerance = .tolerance,
                              info = paste(name, "fitted values", sep = ": "))
    if(.skip[7]) expect_equal(deviance(mod1), deviance(mod2), tolerance = .tolerance,
                              info = paste(name, "deviance", sep = ": "))
    if(.skip[8]) expect_equal(as.numeric(logLik(mod1)), as.numeric(logLik(mod2)), tolerance = .tolerance,
                              info = paste(name, "logLik", sep = ": "))
}


# probit
mode_choice$w <- make_weights(mode_choice)
pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
pbt_glm <- glm(mode ~ cost + ivtime + ovtime, data = mode_choice, family = binomial(link = 'probit'))
check_glm(pbt, pbt_glm, 1E-04, "probit")
pprint_check(pbt, pbt_glm, "probit")

# logit
lgt <- update(pbt, link = 'logit')
lgt_glm <- update(pbt_glm, family = binomial(link = 'logit'))
check_glm(pbt, pbt_glm, 1E-04, "logit")
pprint_check(lgt, lgt_glm, "logit")

# weighted probit
wpbt <- update(pbt, weights = w)
wpbt_glm <- update(pbt_glm, weights = w)
check_glm(wpbt, wpbt_glm, 1E-04, "weighted probit", skip = 8)
pprint_check(wpbt, wpbt_glm, "weighted probit")

# weighted logit
wlgt <- update(lgt, weights = w)
wlgt_glm <- update(lgt_glm, weights = w)
check_glm(wlgt, wlgt_glm, 1E-04, "weighted logit", skip = 8)
pprint_check(wlgt, wlgt_glm, "weighted logit")

# Poisson
trips$w <- make_weights(trips)
eq <- trips ~ workschl + size + dist + smsa + fulltime + distnod + realinc + weekend + car
ps <- poisreg(eq, trips)
ps_glm <- glm(eq, trips, family = poisson)
check_glm(ps, ps_glm, tolerance = 1E-02, name = "Poisson")
pprint_check(ps, ps_glm, "Poisson")

# Negbin2
nb <- poisreg(eq, trips, mixing = "gamma", vlink = "nb2")
nb_glm <- MASS::glm.nb(eq, trips)
check_glm(nb, nb_glm, 1E-02, "Negbin2")
pprint_check(nb, nb_glm, "Negbin2")

# weighted Poisson
wps <- update(ps, weights = w)
wps_glm <- update(ps_glm, weights = w)
check_glm(wps, wps_glm, 1E-02, "weighted Poisson")
pprint_check(wps, wps_glm, "weighted Poisson")

# weighted Negbin2
wnb <- poisreg(eq, trips, mixing = "gamma", vlink = "nb2", weights = w)
wnb_glm <- MASS::glm.nb(eq, trips, weights = w)
check_glm(wnb, wnb_glm, 1E-05, "weighted Negbin2")
pprint_check(wnb, wnb_glm, "weighted Negbin2")

# tobit
charitable$logdon <- with(charitable, log(donation) - log(25))
charitable$w <- runif(nrow(charitable))
eq <- logdon ~ log(income) + religion
tbt_aer <- AER::tobit(eq, data = charitable)
tbt_micsr <- tobit1(eq, data = charitable)
tbt_vgam <- VGAM::vglm(eq, data = charitable, family = VGAM::tobit())
expect_equal(coef(tbt_micsr, subset = "covariates"), coef(tbt_aer), info = "tobit coef")
expect_equal(as.numeric(coef(tbt_micsr, subset = "vcov")), tbt_aer$scale, info = "tobit sigma")
expect_equal(as.numeric(logLik(tbt_aer)), as.numeric(logLik(tbt_micsr)), info = "tobit logLik")
#expect_equal(deviance(tbt_aer), deviance(tbt_micsr), info = "tonit deviance")

tfit <- survreg(Surv(logdon, logdon>0, type='left') ~ log(income) + religion, data=charitable, dist='gaussian')


# weighted tobit

#wtbt_aer <- update(tbt_aer, weights = w)
#wtbt_micsr <- update(tbt_micsr, weights = w)
#expect_equal(coef(wtbt_micsr, subset = "covariates"), coef(wtbt_aer), info = "weighted tobit coef")
#expect_equal(coef(wtbt_micsr, subset = "vcov"), wtbt_aer$scale, info = "weighted tobit sigma")
#expect_equal(as.numeric(logLik(wtbt_aer)), as.numeric(logLik(wtbt_micsr)), info = "weighted tobit logLik")
#expect_equal(deviance(wtbt_aer), deviance(wtbt_micsr), info = "weighted tobit deviance")

# ordreg
ord_micsr <- ordreg(factor(dindx) ~ rhs1 + catchup, fin_reform, link = "logit")
ord_mass <- MASS::polr(factor(dindx) ~ rhs1 + catchup, fin_reform, method = "logistic")
expect_equal(coef(ord_micsr), coef(ord_mass), tolerance = 1E-05, info = "ordreg coef")
expect_equal(as.numeric(logLik(ord_micsr)), as.numeric(logLik(ord_mass)), tolerance = 1E-05, info = "ordreg coef")
expect_equal(as.numeric(deviance(ord_micsr)), as.numeric(deviance(ord_mass)), tolerance = 1E-05, info = "ordreg coef")
expect_equal(fitted(ord_micsr, type = "probabilities"), fitted(ord_mass), tolerance = 1E-05)

# weighted ordreg
fin_reform$w <- runif(nrow(fin_reform))
word_micsr <- update(ord_micsr, weights = w)
word_mass <- update(ord_mass, weights = w)
expect_equal(coef(word_micsr), coef(word_mass), tolerance = 1E-03, info = "ordreg coef")
expect_equal(as.numeric(logLik(word_micsr)), as.numeric(logLik(word_mass)), tolerance = 1E-05, info = "ordreg logLik")
expect_equal(as.numeric(deviance(word_micsr)), as.numeric(deviance(word_mass)), tolerance = 1E-05, info = "ordreg deviance")
expect_equal(fitted(word_micsr, type = "probabilities"), fitted(word_mass), tolerance = 1E-04, info = "ordreg fitted")


## ra <- MASS::polr(mathlevel ~ sat + language + sex + major + mathcourse + physiccourse + chemistcourse, math, method = "probit")
## ra <- micsr::ordreg(mathlevel ~ language + sex + major + mathcourse, math, link = "probit")


## Ey <- function(mu){
##     J <- length(mu)
##     c(- dnorm(mu[1]) / pnorm(mu[1]),
##     (dnorm(mu[1:(J - 1)]) - dnorm(mu[2:J])) / (pnorm(mu[2:J]) -pnorm(mu[1:(J - 1)])),
##     dnorm(mu[J]) / (1 - pnorm(mu[J])))
## }
## Expy <- Ey(ra$zeta)
## y <- model.response(model.frame(ra))
## names(Expy) <- levels(y)
## Expy <- Expy[as.character(y)]
## resid <- Expy - ra$lp

## N <- 1E03
## x <- rnorm(N, sd = 3, mean = 1)
## eps <- rnorm(N)
## alpha <- 1
## beta <- 1
## bx <- alpha + beta * x + eps
## zeta <- c(- Inf, 0, 2, 4, + Inf)
## y <- as.numeric(cut(bx, zeta))
## sim <- tibble(x = x, y = y)
## m <- MASS::polr(factor(y) ~ x, sim, method = "probit")
## ordreg(y ~ x, sim, link = "probit")

## z$w <- runif(nrow(z))
## z$y <- as.numeric(as.factor(z$dindx))
## z$y <- dplyr::case_when(z$y < 5 ~ 1, z$y == 5 ~ 2, z$y == 6 ~ 3, z$y > 6 ~ 4)
## ord_micsr <- ordreg(factor(y) ~ rhs1 + catchup, z, link = "probit")
## ord_mass <- MASS::polr(factor(y) ~ ., z, method = "probit")
## mu <- coef(ord_micsr, subset = "threshold")

## Ey <- function(mu){
##     J <- length(mu)
##     c(- dnorm(mu[1]) / pnorm(mu[1]),
##     (dnorm(mu[1:(J - 1)]) - dnorm(mu[2:J])) / (pnorm(mu[2:J]) -pnorm(mu[1:(J - 1)])),
##     dnorm(mu[J]) / (1 - pnorm(mu[J])))
## }

## Expy <- Ey(mu) 
## names(Expy) <- levels(model.response(model.frame(ord_micsr)))
## Expy <- Expy[as.character(model.response(model.frame(ord_micsr)))]
## res1 <- Expy - ord_micsr$linear.predictors

## res2 <- resids(ord_mass)
## res3 <- presid(ord_mass)

## y <- as.integer(model.response(model.frame(ord_mass)))
## #lp <- ord_micsr$linear.predictor
## lp <- ord_mass$lp
## n <- length(y)
## #zeta <- coef(ord_mass, subset = "threshold")
## zeta <- ord_mass$zeta
## q <- length(zeta)
## pfun <- pnorm
## cumpro <- cbind(0, matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - lp), , q), 1)
## cumpr <- cbind(0, pfun(matrix(zeta, n, q, byrow = TRUE) - lp), 1)
## identical(cumpr, cumpro)
## lo <- cumpr[cbind(seq_len(n), y)]
## hi <- 1 - cumpr[cbind(seq_len(n), y + 1L)]
## res <- lo - hi


## #mu_j - beta ^ x

## U <- matrix(c(- Inf, zeta, + Inf), n, q + 2, byrow = TRUE) - lp
## gres <- (dnorm(U[cbind(seq_len(n), y)]) -
##          dnorm(U[cbind(seq_len(n), y + 1)])) /
##     (pnorm(U[cbind(seq_len(n), y + 1)]) -
##      pnorm(U[cbind(seq_len(n), y)]))


## object <- ord_mass
## n <- length(object$lp)
## q <- length(object$zeta)
## cumpr <- cbind(0, matrix(pfun(matrix(object$zeta, n, q, byrow = TRUE) - 
##                               object$lp), , q), 1)
## y <- as.integer(model.response(object$model))
## lo <- cumpr[cbind(seq_len(n), y)]
## hi <- 1 - cumpr[cbind(seq_len(n), y + 1L)]
## resm <- lo - hi
## res3 <- presid(ord_mass)


## Zeta <- c(-Inf, zeta, +Inf)
## i <- 1:20 ; unname(pnorm(Zeta[y[i]] - lp[i]) + pnorm(Zeta[y[i] + 1] - lp[i]) - 1)

## z <- fin_reform %>% na.omit
## z$w <- runif(nrow(z))
## z$y <- as.numeric(as.factor(z$dindx))
## z$y <- dplyr::case_when(z$y < 5 ~ 1, z$y == 5 ~ 2, z$y > 5 ~ 3)
## m <- MASS::polr(factor(y) ~ gdpg + bop + bank, z, method = "probit")
## l <- ordreg(y ~ gdpg + bop + bank, z, link = "probit")

## load("~/math.rda")
## z <- math
## z$y <- as.numeric(z$mathlevel)
## m <- MASS::polr(factor(y) ~ sat + language + sex + major + mathcourse + physiccourse + chemistcourse, z, method = "probit")
## # checker pourquoi ça plante
## mm <- ordreg(factor(y) ~ sat + language + sex + major + mathcourse + physiccourse + chemistcourse, z, link = "probit")
 

## z <- sim
## y <- z$y
## Zeta <- c(-100, m$zeta, +100)
## n <- length(y)
## Zeta2 <- matrix(Zeta, nrow = n, ncol = length(Zeta), byrow = TRUE)
## Zeta_j <- Zeta2[cbind(1:n, y + 1)]
## Zeta_m1 <- Zeta2[cbind(1:n, y)]
## mu0 <- mean(m$lp)
## mu <- m$lp
## Eeps <- mu0 + (dnorm(Zeta_m1 - mu0) - dnorm(Zeta_j - mu0)) / (pnorm(Zeta_j - mu0) - pnorm(Zeta_m1 - mu0))
## Eepsx <- mu + (dnorm(Zeta_m1 - mu) -  dnorm(Zeta_j - mu))  / (pnorm(Zeta_j - mu)  - pnorm(Zeta_m1 - mu))
## res <- Eeps - mu
## resx <- Eepsx - mu




## # recommandations muscler la methodo et recommandations politique éco




## N <- 1E03
## x <- rnorm(N, sd = 3, mean = 1)
## eps <- rnorm(N)
## alpha <- 1
## beta <- 1
## bx <- alpha + beta * x + eps
## zeta <- c(- Inf, 0, 2, 4, + Inf)
## y <- as.numeric(cut(bx, zeta))
## sim <- tibble(x = x, y = y)
## m <- MASS::polr(factor(y) ~ x, sim, method = "probit")
## o <- ordreg(y ~ x, sim, link = "probit")
## zeta <- coef(o, subset = "threshold")
## q <- length(zeta)
## U <- matrix(c(- Inf, zeta), n, q + 2, byrow = TRUE) -
##     o$linear.predictor
## gres <- (dnorm(U[cbind(seq_len(n), y)]) -
##          dnorm(U[cbind(seq_len(n), y + 1)])) /
##     (pnorm(U[cbind(seq_len(n), y + 1)]) -
##      pnorm(U[cbind(seq_len(n), y)]))



## U <- matrix(c(- Inf, zeta, + Inf), n, q + 2, byrow = TRUE) - lp
## gres <- (dnorm(U[cbind(seq_len(n), y)]) -
##          dnorm(U[cbind(seq_len(n), y + 1)])) /
##     (pnorm(U[cbind(seq_len(n), y + 1)]) -
##      pnorm(U[cbind(seq_len(n), y)]))
