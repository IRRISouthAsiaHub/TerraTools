## -----------------------------------------------------------------------------
#| warning: false
library(micsr)
bank <- binomreg(federiv ~ eqrat + optval +  mktbk +
                      perfor + dealdum | . - eqrat - optval +
                      no_emp + no_subs + no_off,
                  data = federiv, method = "ml")


## -----------------------------------------------------------------------------
names(bank$coefficients)


## -----------------------------------------------------------------------------
bank$npar


## -----------------------------------------------------------------------------
coef(bank)


## -----------------------------------------------------------------------------
coef(bank, subset = c("chol", "resid"))


## -----------------------------------------------------------------------------
#| eval: false
# coef(bank, subset = NULL)


## -----------------------------------------------------------------------------
coef(bank, coef = c("mktbk", "perfor"))


## -----------------------------------------------------------------------------
coef(bank, grep = "eqrat|perfor")


## -----------------------------------------------------------------------------
coef(bank, grep = "eqrat|perfor", subset = c("instruments", "chol"))


## -----------------------------------------------------------------------------
coef(bank, subset = "instruments", grep = "Intercept", invert = TRUE)


## -----------------------------------------------------------------------------
vcov(bank, subset = "resid")


## -----------------------------------------------------------------------------
vcov(bank, subset = "resid")
vcov(bank, subset = "resid", vcov = "hessian")
vcov(bank, subset = "resid", vcov = "opg")
vcov(bank, subset = "resid", vcov = "hc")


## -----------------------------------------------------------------------------
stder(bank, subset = "resid")
stder(bank, vcov = "hc", subset = "resid")


## -----------------------------------------------------------------------------
charitable |> dummy( education, religion) |> head(2)


## -----------------------------------------------------------------------------
charitable |> dummy(religion, keep = TRUE) |> head(2)
charitable |> dummy(religion, ref = TRUE) |> head(2)


## -----------------------------------------------------------------------------
t.test(donation ~ married, charitable)


## -----------------------------------------------------------------------------
#| collapse: true
t.test(donation ~ married, charitable) |> gaze()


## -----------------------------------------------------------------------------
tobit1(log(donation / 25) ~ married, charitable) |> gaze()

