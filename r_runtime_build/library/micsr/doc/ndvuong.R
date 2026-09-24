## -----------------------------------------------------------------------------
#| collapse: true
library(micsr)
Vuong <- vuong_sim(N = 100, R = 1000, Kf = 15, Kg = 1, a = 0.5)
head(Vuong)
mean(Vuong)
mean(abs(Vuong) > 1.96)


## -----------------------------------------------------------------------------
#| label: fig-vuong
#| fig-cap: "Empirical distribution if the Vuong statistic"
#| echo: false
#| message: false
if (requireNamespace("ggplot2")){
    library(ggplot2)
    ggplot(data = data.frame(Vuong = Vuong)) + geom_density(aes(x = Vuong)) +
        geom_function(fun = dnorm, linetype = "dotted")
}


## -----------------------------------------------------------------------------
test <- ndvuong(turnout$group, turnout$intens, size = 0.05, pval = FALSE)
test


## -----------------------------------------------------------------------------
#| collapse: true
ndvuong(turnout$group, turnout$intens, nd = FALSE) |> gaze()


## -----------------------------------------------------------------------------
test <- ndvuong(turnout$group, turnout$intens, pval = TRUE)
test


## -----------------------------------------------------------------------------
#| message: false
if (requireNamespace("mlogit")){
    library(mlogit)
    data("ModeCanada", package = "mlogit")
    MC <- dfidx(ModeCanada, subset = noalt == 4)
}


## -----------------------------------------------------------------------------
if (requireNamespace("mlogit")){
    ml <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
                 reflevel = 'car', alt.subset = c("car", "train", "air"))
}


## -----------------------------------------------------------------------------
if (requireNamespace("mlogit")){
    hl <- mlogit(choice ~ freq + cost + ivt + ovt | urban + income, MC, 
                 reflevel = 'car', alt.subset = c("car", "train", "air"),
                 heterosc = TRUE)
    coef(summary(hl))
}


## -----------------------------------------------------------------------------
if (requireNamespace("lmtest")){
    lmtest::waldtest(hl, heterosc = FALSE) |> gaze()
    scoretest(ml, heterosc = TRUE) |> gaze()
    lmtest::lrtest(hl, ml) |> gaze()
}


## -----------------------------------------------------------------------------
ndvuong(hl, ml, nested = TRUE)


## -----------------------------------------------------------------------------
if (requireNamespace("mlogit")){
    library(mlogit)
    data("RiskyTransport", package = "mlogit")
    RT <- dfidx(RiskyTransport, idx = c(id = "chid", "mode"),
                choice = "choice")
}


## -----------------------------------------------------------------------------
if (requireNamespace("mlogit")){
    ml <- mlogit(choice ~ cost, data = RT)
    hl <- mlogit(choice ~ cost , data = RT, heterosc = TRUE)
    nl <- mlogit(formula = choice ~ cost, data = RT,
                 nests = list(fast = c("Helicopter", "Hovercraft"),
                              slow = c("WaterTaxi", "Ferry")),
                 un.nest.el = TRUE)
    xl <- mlogit(choice ~ cost, data = RT, rpar = c(cost = "zbt"))
}


## -----------------------------------------------------------------------------
#| collapse: true
ndvuong(nl, hl, vartest = TRUE) |> gaze()


## -----------------------------------------------------------------------------
ndvuong(hl, nl)


## -----------------------------------------------------------------------------
ndvuong(xl, ml)

