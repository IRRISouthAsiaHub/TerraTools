# test.partykit.R: test partykit and evtree packages

source("test.prolog.R")
library(plotmo)
library(partykit)
data("BostonHousing", package = "mlbench")
data("PimaIndiansDiabetes", package = "mlbench")

# lmtree

boston <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))
set.seed(2018)
lmtree.boston1 <- lmtree(medv ~ log(lstat) + rm^2 |
  crim + ptratio + tax + dis + rad + chas,
  data = boston, minsize = 40)

boston2 <- boston
boston2$log.lstat <- log(boston2$lstat)
boston2$lstat <- NULL
boston2$rm.squared <- boston2$rm^2
boston2$rm <- NULL
set.seed(2018)
lmtree.boston2 <- lmtree(medv ~ log.lstat + rm.squared |
  crim + ptratio + tax + dis + rad + chas,
  data = boston2, minsize = 40)

plot(lmtree.boston1)
plot(lmtree.boston2)

plotmo(lmtree.boston1, SHOWCALL=TRUE)
plotmo(lmtree.boston2, trace=2, SHOWCALL=TRUE)
plotmo(lmtree.boston2, trace=1, all1=TRUE, degree2=c("ptratio", "log.lstat"), SHOWCALL=TRUE)
plotmo(lmtree.boston2, all1=TRUE, all2=TRUE, SHOWCALL=TRUE)

# TODO gives warnings because of because of price/citations in formula
# data("Journals", package = "AER")
# Journals <- transform(Journals,
#   age = 2000 - foundingyear,
#   chars = charpp * pages)
# j_tree <- lmtree(log(subs) ~ log(price/citations) | price + citations +
#   age + chars + society, data = Journals, minsize = 10)
# plotmo(j_tree, SHOWCALL=TRUE)

# Works, but commented out to save testing time:
# data("TeachingRatings", package = "AER")
# tr_tree <- lmtree(eval ~ beauty | age + gender + division,
#    data = TeachingRatings, weights = students, subset = credits == "more",
#    caseweights = FALSE)
# plot(tr_tree)
# plotmo(tr_tree, all1=TRUE, all2=TRUE, SHOWCALL=TRUE)

# glmtree

glmtree1 <- glmtree(diabetes ~ glucose | mass + age,
                    data = PimaIndiansDiabetes, family = binomial)
plot(glmtree1)
plotmo(glmtree1, SHOWCALL=TRUE)
plotmo(glmtree1, all2=TRUE, SHOWCALL=TRUE)

# mob

pima <- PimaIndiansDiabetes[1:50,] # small set of data for fast test

logit1 <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    # note that a complicated formula is necessary
    formula <- as.formula(paste("y ~ ", paste(colnames(x)[-1], collapse="+"))) # -1 drops intercept
    glm(formula=formula, data=as.data.frame(x), family=binomial, start=start, ...)
}
mob1 <- mob(diabetes ~ glucose | mass + age,
            data = PimaIndiansDiabetes, fit = logit1)
plot(mob1)
plotmo(mob1, trace=1, SHOWCALL=TRUE)
plotmo(mob1, pmethod="partdep", degree1=0,
       degree2=c("glucose", "mass"), persp.ticktype="detailed", SHOWCALL=TRUE)
plotmo(mob1, all1=TRUE, all2=TRUE, SHOWCALL=TRUE)

logit2 <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    glm(y ~ 0 + x, family = binomial, start = start, ...)
}
mob2 <- mob(diabetes ~ glucose | mass, data = pima, fit = logit2)
expect.err(try(plotmo(mob2)), "The formula in the mob fit function is not supported by plotmo")

logit3 <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    glm(y ~ 0+x  , family = binomial, start = start, ...)
}
mob3 <- mob(diabetes ~ glucose | age, data = pima, fit = logit3)
expect.err(try(plotmo(mob3)), "The formula in the mob fit function is not supported by plotmo")

logit4 <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    glm(y ~ x - 1, family = binomial, start = start, ...)
}
mob4 <- mob(diabetes ~ glucose | age, data = pima, fit = logit4)
expect.err(try(plotmo(mob4)), "The formula in the mob fit function is not supported by plotmo")

logit5 <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    glm(y~x-1 , family = binomial, start = start, ...)
}
mob5 <- mob(diabetes ~ glucose | age, data = pima, fit = logit5)
expect.err(try(plotmo(mob5)), "The formula in the mob fit function is not supported by plotmo")

logit6 <- function (y, x, start = NULL, weights = NULL, offset = NULL, ...)
{
    glm(as.formula(paste("y ~ ", paste(colnames(x)[-1], collapse="+"))),
        data=data.frame(x), family = binomial, start = start, ...)
}
mob6 <- mob(diabetes ~ glucose | mass + age, data = pima, fit = logit6)
plot(mob6) # tree is just a root (no branches)
plotmo(mob6)

library(rpart.plot)
rpart.Kyphosis <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
plotmo(rpart.Kyphosis, SHOWCALL=TRUE)
party.Kyphosis <- as.party(rpart.Kyphosis)
expect.err(try(plotmo(party.Kyphosis)), "cannot get the original model predictors")

library(evtree)
## regression
set.seed(1090)
airq <- subset(airquality, !is.na(Ozone) & complete.cases(airquality))
ev_air <- evtree(Ozone ~ ., data = airq)
# plot(ev_air)
plotmo(ev_air, SHOWCALL=TRUE)
## classification
ev_iris <- evtree(Species ~ .,data = iris)
# plot(ev_iris)
plotmo(ev_iris, SHOWCALL=TRUE)
plotmo(ev_iris, type="prob", nresponse="versicolor", pmethod="apartdep", SHOWCALL=TRUE)
plotres(ev_iris, type="prob", nresponse="setosa", SHOWCALL=TRUE)

# cforest

cforest1 <- cforest(dist ~ speed, data = cars)
plotmo(cforest1, trace=1, SHOWCALL=TRUE)
plotres(cforest1, trace=1, SHOWCALL=TRUE)

data("mammoexp", package = "TH.data")
cforest2 <- cforest(ME ~ PB + SYMPT, data = mammoexp, ntree = 5)
plotmo(cforest2, trace=1, SHOWCALL=TRUE, pmethod="apartdep")
plotres(cforest2)

source("test.epilog.R")
