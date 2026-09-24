## -----------------------------------------------------------------------------
#| warning: false
library(micsr)
bank_msq <- ivldv(federiv ~ eqrat + optval +  mktbk +
                      perfor + dealdum | . - eqrat - optval +
                      no_emp + no_subs + no_off,
                  data = federiv, method = "ml")


## -----------------------------------------------------------------------------
bank_msq$npar


## -----------------------------------------------------------------------------
#| collapse: true
select_coef(bank_msq)
select_coef(bank_msq, subset = c("resid", "chol"))


## -----------------------------------------------------------------------------
#| collapse: true
coef(bank_msq)
coef(bank_msq, subset = c("resid", "chol"))
vcov(bank_msq, subset = c("resid"))


## -----------------------------------------------------------------------------
#| collapse: true
vcov(bank_msq, subset = "all", grep = "Intercept")


## -----------------------------------------------------------------------------
#| collapse: true
npar(bank_msq)
npar(bank_msq, subset = c("resid", "chol"))


## -----------------------------------------------------------------------------
vcov(bank_msq, subset = "resid", vcov = "opg")


## -----------------------------------------------------------------------------
#| collapse: true
sqrt(diag(vcov(bank_msq, subset = "resid")))
stder(bank_msq, subset = "resid")


## -----------------------------------------------------------------------------
pbt <- binomreg(mode ~ cost + ivtime + ovtime,
                data = mode_choice, link = 'probit')
pbt$logLik


## -----------------------------------------------------------------------------
#| collapse: true
logLik(pbt)
logLik(pbt, type = "model")
logLik(pbt, type = "saturated")
logLik(pbt, type = "null")


## -----------------------------------------------------------------------------
#| collapse: true
AIC(pbt)
BIC(pbt)
AIC(pbt, type = "null")
AIC(pbt, k = 5)


## -----------------------------------------------------------------------------
#| collapse: true
deviance(pbt)
deviance(pbt, type = "null")


## -----------------------------------------------------------------------------
#| collapse: true
pbt$test


## -----------------------------------------------------------------------------
summary(bank_msq, subset = c("chol", "resid"), vcov = "opg")

