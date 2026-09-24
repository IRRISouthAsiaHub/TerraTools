## -----------------------------------------------------------------------------
#| echo: true
#| message: false
library(micsr)
random_group |> head()


## -----------------------------------------------------------------------------
random_group <- transform(random_group, treated = ifelse(group < 0, 1, 0))


## -----------------------------------------------------------------------------
mod0 <- lm(log(wage) ~ treated + age + child, data = random_group)


## -----------------------------------------------------------------------------
names(mod0)


## -----------------------------------------------------------------------------
#| collapse: true
mod0$rank
length(mod0$coefficients)


## -----------------------------------------------------------------------------
mod0$terms


## -----------------------------------------------------------------------------
#| collapse: true
mod0$assign


## -----------------------------------------------------------------------------
#| collapse: true
coef(mod0)
df.residual(mod0)
resid(mod0) |> head()
fitted(mod0) |> head()


## -----------------------------------------------------------------------------
#| results: false
terms(mod0)


## -----------------------------------------------------------------------------
#| collapse: true
vcov(mod0)
sigma(mod0)
nobs(mod0)


## -----------------------------------------------------------------------------
head(mod0$model, 3)


## -----------------------------------------------------------------------------
#| results: false
model.frame(mod0)


## -----------------------------------------------------------------------------
#| collapse: true
names(attributes(mod0$model))


## -----------------------------------------------------------------------------
model.matrix(terms(mod0), model.frame(mod0)) |> head(3)


## -----------------------------------------------------------------------------
#| results: hide
#| collapse: true
model.matrix(mod0)


## -----------------------------------------------------------------------------
dim(model.matrix(mod0))


## -----------------------------------------------------------------------------
#| collapse: true
mod0 |> model.frame() |> model.response() |> head()


## -----------------------------------------------------------------------------
#| collapse: true
class(mod0$call)
mod0$call


## -----------------------------------------------------------------------------
mod1 <- lm(log(wage) ~ treated + age + child + female + single + migrant +
               temp, data = random_group)


## -----------------------------------------------------------------------------
mod1 <- mod0 |> update(log(wage) ~ treated + age + child +
                            female + single + migrant + temp)


## -----------------------------------------------------------------------------
mod1 <- update(mod0, . ~ . + female + single + migrant + temp)


## -----------------------------------------------------------------------------
#| eval: false
# lm(log(wage) ~ treated + age + child + female + single + migrant + temp,
#            data = random_group, subset = age >= 20)


## -----------------------------------------------------------------------------
#| eval: false
# update(mod0, . ~ . + female + single + migrant + temp, subset = age >= 20)


## -----------------------------------------------------------------------------
#| collapse: true
options()$na.action


## -----------------------------------------------------------------------------
#| results: false
mod1 |> update(na.action = "na.fail")


## -----------------------------------------------------------------------------
#| eval: false
# mod1 |> update(. ~ . + ten, na.action = "na.fail")


## -----------------------------------------------------------------------------
mod2 <- mod1 |> update(. ~ . + ten)
names(mod2)


## -----------------------------------------------------------------------------
unname(mod2$na.action)


## -----------------------------------------------------------------------------
#| collapse: true
mod2$model |> attributes() |> names()


## -----------------------------------------------------------------------------
#| collapse: true
levels(random_group$fsize)
levels(random_group$edu)


## -----------------------------------------------------------------------------
mod3 <- lm(log(wage) ~ treated + poly(age, 2) + child + fsize + edu + 
               female + single + migrant + temp + ten,
            data = random_group)


## -----------------------------------------------------------------------------
mod3$assign


## -----------------------------------------------------------------------------
names(coef(mod3))


## -----------------------------------------------------------------------------
mod3$xlevels


## -----------------------------------------------------------------------------
mod3$contrasts


## -----------------------------------------------------------------------------
mod4 <- update(mod3, contrasts = list(edu = "contr.sum",
                                      fsize = "contr.helmert"))
mod4$contrasts


## -----------------------------------------------------------------------------
head(mod4$model, 3)


## -----------------------------------------------------------------------------
head(model.matrix(mod3), 3)
head(model.matrix(mod4), 3)


## -----------------------------------------------------------------------------
mod4 <- update(mod3, weights = samplew)


## -----------------------------------------------------------------------------
model.frame(mod4) |> head(3)


## -----------------------------------------------------------------------------
model.weights(model.frame(mod4)) |> head(3)


## -----------------------------------------------------------------------------
random_group <- transform(random_group, v = 3/4 + 1/4 * as.numeric(edu))
mod5 <- update(mod4, . ~ . + offset(v) - edu)
mod6 <- update(mod4, . ~ . - edu, offset = v)


## -----------------------------------------------------------------------------
model.offset(model.frame(mod5)) |> head(3)


## -----------------------------------------------------------------------------
mylm <- function (formula, data, subset, weights, na.action, method = "qr", 
                  model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                  singular.ok = TRUE, contrasts = NULL, offset, ...){
    mf <- match.call()
    mf
}


## -----------------------------------------------------------------------------
mf <- mylm(log(wage) ~ treated + poly(age, 2) + child + fsize + female +
               single + migrant + temp + ten,
           random_group,
           subset = age > 20,
           weights = samplew,
           contrasts = list(fsize = "contr.helmert"),
           offset = v
           )
mf
cl <- mf


## -----------------------------------------------------------------------------
names(mf)
m <- match(c("formula", "data", "subset", "weights", "na.action", 
             "offset"), names(mf), 0L)
m


## -----------------------------------------------------------------------------
mf <- mf[c(1L, m)]
mf


## -----------------------------------------------------------------------------
mf[[1L]] <- quote(stats::model.frame)
mf


## -----------------------------------------------------------------------------
mf <- eval(mf)
head(mf, 3)
names(attributes(mf))


## -----------------------------------------------------------------------------
mt <- attr(mf, "terms")                                            
y <- model.response(mf)                                 
w <- model.weights(mf)                                  
offset <- model.offset(mf)                                         
x <- model.matrix(mt, mf, list(fsize = "contr.helmert"))


## -----------------------------------------------------------------------------
z <- lm.wfit(x, y, w, offset = offset)
names(z)


## -----------------------------------------------------------------------------
class(z) <- "lm"


## -----------------------------------------------------------------------------
z$na.action <- attr(mf, "na.action")
z$offset <- offset
z$contrasts <- attr(x, "contrasts")
z$xlevels <- .getXlevels(mt, mf)
z$call <- cl
z$terms <- mt
z$model <- mf


## -----------------------------------------------------------------------------
z


## -----------------------------------------------------------------------------
mylm <- function (formula, data, subset, weights, na.action, method = "qr", 
                  model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                  singular.ok = TRUE, contrasts = NULL, offset, ...){
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")                                            
    y <- model.response(mf, "numeric")                                 
    w <- as.vector(model.weights(mf))                                  
    offset <- model.offset(mf)                                         
    x <- model.matrix(mt, mf, contrasts)
    z <- lm.wfit(x, y, w, offset = offset)
    class(z) <- "lm"
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$model <- mf
    z
}


## -----------------------------------------------------------------------------
mod6 <- mylm(log(wage) ~ treated + poly(age, 2) + child + fsize + female +
                 single + migrant + temp + ten, 
             random_group,
             subset = age > 20,
             weights = samplew,
             contrasts = list(fsize = "contr.helmert"),
             offset = v)

