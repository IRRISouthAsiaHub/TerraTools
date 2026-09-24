# test.plotmo3.R: extra tests for plotmo version 3 and higher

source("test.prolog.R")
library(earth)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur

# check check.numeric.scalar

xtest <- NA
expect.err(try(plotmo:::check.numeric.scalar(xtest)), "'xtest' is NA")
xtest <- NULL
expect.err(try(plotmo:::check.numeric.scalar(xtest)), "'xtest' is NULL")
expect.err(try(plotmo:::check.numeric.scalar(NA)), "argument is NA")
expect.err(try(plotmo:::check.numeric.scalar(NULL)), "argument is NULL")
expect.err(try(plotmo:::check.numeric.scalar(try)), "'try' must be numeric (whereas its current class is \"function\")")
expect.err(try(plotmo:::check.numeric.scalar('try')), "\"try\" must be numeric (whereas its current class is \"character\")")
expect.err(try(plotmo:::check.numeric.scalar(NULL)), "argument is NULL")
expect.err(try(plotmo:::check.numeric.scalar(1234, min=2, max=3)), "argument=1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar(0.1234, min=2, max=3)), "argument=0.1234 but it should be between 2 and 3")

expect.err(try(plotmo:::check.numeric.scalar(.1234, min=2, max=3)), "argument=0.1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar(+1234, min=2, max=3)), "argument=1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar(-1234, min=2, max=3)), "argument=-1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar(+.1234, min=2, max=3)), "argument=0.1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar(-.1234, min=2, max=3)), "argument=-0.1234 but it should be between 2 and 3")
expect.err(try(plotmo:::check.numeric.scalar("", min=0, max=3)), "\"\" must be numeric (whereas its current class is \"character\"")

x.numeric.scalar <- 1234
expect.err(try(plotmo:::check.numeric.scalar(x.numeric.scalar, min=0, max=3)), "x.numeric.scalar=1234 but it should be between 0 and 3")
stopifnot(identical(plotmo:::check.numeric.scalar(x.numeric.scalar, min=2, max=1235), 1234))
stopifnot(identical(plotmo:::check.numeric.scalar(1234, min=2, max=1235), 1234))

# check check.integer.scalar

xtest <- NA
expect.err(try(plotmo:::check.integer.scalar(xtest)), "'xtest' is NA")
xtest <- NULL
expect.err(try(plotmo:::check.integer.scalar(xtest)), "'xtest' is NULL")
expect.err(try(plotmo:::check.integer.scalar(NA)), "argument is NA")
expect.err(try(plotmo:::check.integer.scalar(NA, null.ok=TRUE)), "argument is NA")
expect.err(try(plotmo:::check.integer.scalar(NULL)), "argument is NULL")
expect.err(try(plotmo:::check.integer.scalar(xtest, na.ok=TRUE)), "'xtest' is NULL")
expect.err(try(plotmo:::check.integer.scalar("xyz", na.ok=TRUE)), "\"xyz\" is a string but it should be an integer, or NA, or TRUE or FALSE")
expect.err(try(plotmo:::check.integer.scalar("TRUE", na.ok=TRUE)), "\"TRUE\" is a string but it should be an integer, or NA, or TRUE or FALSE")
stopifnot(identical(plotmo:::check.integer.scalar(TRUE), TRUE))
stopifnot(identical(plotmo:::check.integer.scalar(NA, na.ok=TRUE), NA))
x.integer.scalar <- 1234L
expect.err(try(plotmo:::check.integer.scalar(x.integer.scalar, min=0, max=3)), "x.integer.scalar=1234 but it should be between 0 and 3")
stopifnot(identical(plotmo:::check.integer.scalar(x.integer.scalar, min=2, max=1235), 1234L))
stopifnot(identical(plotmo:::check.integer.scalar(1234, min=2, max=1235), 1234))
stopifnot(identical(plotmo:::check.integer.scalar(x.integer.scalar, min=2, max=1235), 1234L))
stopifnot(identical(plotmo:::check.integer.scalar(1234, min=2, max=1235), 1234))
xtest <- 1.234
expect.err(try(plotmo:::check.integer.scalar(xtest, min=0, max=3)), "xtest=1.234 but it should be an integer, or TRUE or FALSE")

# check check.vec
xtest <- "x"
expect.err(try(plotmo:::check.vec(xtest, "xtest", na.ok=TRUE)), "'xtest' is not numeric")
xtest <- as.double(NA)
print(plotmo:::check.vec(xtest, "xtest", na.ok=TRUE))
xtest <- as.double(1:3)
print(plotmo:::check.vec(xtest, "xtest", na.ok=TRUE))
xtest <- c(1,2,3,1/0,5,6,7)
expect.err(try(plotmo:::check.vec(xtest, "xtest", na.ok=TRUE)), "non-finite value in xtest")
xtest <- c(1,2,3,NA,5,6,7)
expect.err(try(plotmo:::check.vec(xtest, "xtest")), "NA in xtest")
xtest <- c(1,2,3)
expect.err(try(plotmo:::check.vec(xtest, "xtest", expected.len=2)), "'xtest' has the wrong length 3, expected 2")
print(plotmo:::check.vec(c(TRUE, FALSE), "c(TRUE, FALSE)"))

plotmo1 <- function(object, ..., trace=0, SHOWCALL=TRUE, caption=NULL) {
    if(is.null(caption))
        caption <- paste(deparse(substitute(object)), collapse=" ")
    call <- match.call(expand.dots=TRUE)
    call <- strip.space(paste(deparse(substitute(call)), collapse=" "))
    printf("%s\n", call)
    plotmo(object, trace=trace, SHOWCALL=SHOWCALL, caption=caption, ...)
}
plotres1 <- function(object, ..., trace=0, SHOWCALL=TRUE, caption=NULL) {
    if(is.null(caption))
        caption <- paste(deparse(substitute(object)), collapse=" ")
    call <- match.call(expand.dots=TRUE)
    call <- strip.space(paste(deparse(substitute(call)), collapse=" "))
    printf("%s\n", call)
    plotres(object, trace=trace, SHOWCALL=SHOWCALL, caption=caption, ...)
}
# basic tests of plotmo on abbreviated titanic data

get.tita <- function()
{
    tita <- etitanic
    pclass <- as.character(tita$pclass)
    # change the order of the factors so not alphabetical
    pclass[pclass == "1st"] <- "first"
    pclass[pclass == "2nd"] <- "class2"
    pclass[pclass == "3rd"] <- "classthird"
    tita$pclass <- factor(pclass, levels=c("class2", "classthird", "first"))
    # log age is so we have a continuous predictor even when model is age~.
    set.seed(2015)
    tita$logage <- log(tita$age) + rnorm(nrow(tita))
    tita$parch <- NULL
    # by=12 gives us a small fast model with an additive and a interaction term
    tita[seq(1, nrow(etitanic), by=12), ]
}
tita <- get.tita()

mod.lm.age <- lm(age~., data=tita)
plotmo1(mod.lm.age)
plotmo1(mod.lm.age, level=.95)
plotmo1(mod.lm.age, level=.95, col.resp=3)

sexn <- as.numeric(tita$sex)
mod.lm.sexn <- lm(sexn~.-sex, data=tita)
plotmo1(mod.lm.sexn)
plotmo1(mod.lm.sexn, level=.95)

set.seed(2020)
mod.earth.age <- earth(age~., data=tita, degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(mod.earth.age)
plotmo1(mod.earth.age, level=.9, degree2=0)

# tita[,4] is age
set.seed(2020)
mod.earth.tita.age <- earth(tita[,-4], tita[,4], degree=2, nfold=3, ncross=3, trace=.5, varmod.method="lm")
cat("\nsummary(mod.earth.tita.age)\n")
print(summary(mod.earth.tita.age))
plotmo1(mod.earth.tita.age)
plotmo1(mod.earth.tita.age, level=.9, degree2=0)

set.seed(2020)
a.earth.sex <- earth(sex~., data=tita, degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(a.earth.sex)
plotmo1(a.earth.sex, level=.9)
plotmo1(a.earth.sex, type="class")
expect.err(try(plotmo1(a.earth.sex, level=.9, degree2=0, type="class")), "predicted values are strings")

# tita[,3] is sex
set.seed(2020)
mod.earth.tita <- earth(tita[,-3], tita[,3], degree=2, nfold=3, ncross=3, varmod.method="lm")
plotmo1(mod.earth.tita)
plotmo1(mod.earth.tita, level=.9, degree2=0)
plotmo1(mod.earth.tita, type="class")
expect.err(try(plotmo1(mod.earth.tita, level=.9, degree2=0, type="class")), "predicted values are strings")

set.seed(2020)
mod.earth.sex <- earth(sex~., data=tita, degree=2, nfold=3, ncross=3, varmod.method="earth", glm=list(family=binomial))
plotmo1(mod.earth.sex)
plotmo1(mod.earth.sex, type="link")
plotmo1(mod.earth.sex, type="class")
plotmo1(mod.earth.sex, level=.9, type="earth")

# tita[,3] is sex
set.seed(2020)
mod.earth.tita <- earth(tita[,-3], tita[,3], degree=2, nfold=3, ncross=3, varmod.method="earth", glm=list(family=binomial))
plotmo1(mod.earth.tita)
plotmo1(mod.earth.tita, type="link")
plotmo1(mod.earth.tita, type="class")
plotmo1(mod.earth.tita, level=.9, type="earth")

# check factor handling when factors are not ordered alphabetically
tita.orgpclass <- etitanic[seq(1, nrow(etitanic), by=12), ]
tita  <- get.tita()
tita$logage <- NULL
tita.orgpclass$parch <- NULL
stopifnot(names(tita.orgpclass) == names(tita))
a.tita.orgpclass <- earth(pclass~., degree=2, data=tita.orgpclass)
a.tita           <- earth(pclass~., degree=2, data=tita)
options(warn=2) # treat warnings as errors
expect.err(try(plotmo(a.tita)), "Defaulting to nresponse=1, see above messages")
options(warn=1)
# following two graphs should be identical
plotmo1(a.tita.orgpclass, nresponse="1st",   all1=T, col.resp=3, type2="im")
plotmo1(a.tita,           nresponse="first", all1=T, col.resp=3, type2="im")
# following two graphs should be identical
plotmo1(a.tita.orgpclass, nresponse="2nd",    all1=T)
plotmo1(a.tita,           nresponse="class2", all1=T)

tita  <- get.tita()
mod.earth.pclass <- earth(pclass~., data=tita, degree=2)
options(warn=2) # treat warnings as errors
expect.err(try(plotmo1(mod.earth.pclass)), "Defaulting to nresponse=1, see above messages")
options(warn=1)
plotmo1(mod.earth.pclass, nresponse="fi")
plotmo1(mod.earth.pclass, nresponse="first")
plotmo1(mod.earth.pclass, nresponse=3)
plotmo1(mod.earth.pclass, type="class")
plotmo1(mod.earth.pclass, nresponse=1,
       type="class", grid.levels=list(sex="fem"),
       smooth.col="indianred", smooth.lwd=2,
       pt.col=as.numeric(tita$pclass)+1,
       pt.pch=1)

# tita[,1] is pclass
mod.earth.tita <- earth(tita[,-1], tita[,1], degree=2)
options(warn=2) # treat warnings as errors
expect.err(try(plotmo1(mod.earth.tita)), "Defaulting to nresponse=1, see above messages")
options(warn=1)
plotmo1(mod.earth.tita, nresponse="first")
plotmo1(mod.earth.tita, type="class")

mod.earth.pclass2 <- earth(pclass~., data=tita, degree=2, glm=list(family=binomial))
# expect.err(try(plotmo1(mod.earth.pclass2)), "nresponse is not specified")
plotmo1(mod.earth.pclass2, nresponse=3)
plotmo1(mod.earth.pclass2, type="link", nresponse=3)
plotmo1(mod.earth.pclass2, type="class")

# tita[,1] is pclass
mod.earth.tita <- earth(tita[,-1], tita[,1], degree=2, glm=list(family=binomial))
plotmo1(mod.earth.tita, nresponse=3)
plotmo1(mod.earth.tita, type="link", nresponse=3)
plotmo1(mod.earth.tita, type="class")

# plotmo vignette examples

# use a small set of variables for illustration
printf("library(earth)\n")
library(earth) # for ozone1 data
data(ozone1)
oz <- ozone1[, c("O3", "humidity", "temp", "ibt")]

lm.model.vignette <- lm(O3 ~ humidity + temp*ibt, data=oz) # linear model
plotmo1(lm.model.vignette, pt.col="gray", nrug=-1)
plotmo1(lm.model.vignette, level=.9)

printf("library(mda)\n")
library(mda)
mars.model.vignette1 <- mars(oz[,-1], oz[,1], degree=2)
plotmo1(mars.model.vignette1)
plotres1(mars.model.vignette1)
mars.model.vignette2 <- mars(oz[,-1,drop=FALSE], oz[,1,drop=FALSE], degree=2)
plotmo1(mars.model.vignette2)
# TODO causes Error in lm.fit(object$x, y, singular.ok = FALSE) : (list) object cannot be coerced to type 'double'
#      although still works
#      the error is mars.to.earth try(hatvalues.lm.fit(lm.fit(object$x, y, singular.ok=FALSE)))
plotres1(mars.model.vignette2, trace=1)

printf("library(rpart)\n")
library(rpart)                                          # rpart
rpart.model.vignette <- rpart(O3 ~ ., data=oz)
plotmo1(rpart.model.vignette, all2=TRUE)
expect.err(try(plotmo1(rpart.model.vignette, level=.9)), "the level argument is not supported for \"rpart\" objects")

# commented out because is slow and already tested in test.non.earth.R
# printf("library(randomForest)\n")
# library(randomForest)                                   # randomForest
# rf.model.vignette <- randomForest(O3~., data=oz)
# plotmo1(rf.model.vignette)
# partialPlot(rf.model.vignette, oz, temp) # compare to partial-dependence plot

printf("library(gbm)\n")
library(gbm)                                            # gbm
set.seed(2016)
gbm.model.vignette <- gbm(O3~., data=oz, dist="gaussian", inter=2, n.trees=100)
# commented out following because they always take the whole page
# plot(gbm.model.vignette, i.var=2) # compare to partial-dependence plots
# plot(gbm.model.vignette, i.var=c(2,3))
set.seed(2016)
plotmo1(gbm.model.vignette, caption="gbm.model.vignette")

# commented out because is slow and already tested elsewhere
# printf("library(mgcv)\n")
# library(mgcv)                                           # gam
# gam.model.vignette <- gam(O3 ~ s(humidity)+s(temp)+s(ibt)+s(temp,ibt), data=oz)
# plotmo1(gam.model.vignette, level=.95, all2=TRUE)

printf("library(nnet)\n")
library(nnet)                                           # nnet
set.seed(4)
nnet.model.vignette <- nnet(O3~., data=scale(oz), size=2, decay=0.01, trace=FALSE)
plotmo1(nnet.model.vignette, type="raw", all2=T)

printf("library(MASS)\n")
library(MASS)                                           # qda
lcush <- data.frame(Type=as.numeric(Cushings$Type),log(Cushings[,1:2]))
lcush <- lcush[1:21,]
qda.model.vignette <- qda(Type~., data=lcush)
plotmo1(qda.model.vignette, type="class", all2=TRUE,
       type2="contour", ngrid2=100, contour.nlevels=2, contour.drawlabels=FALSE,
       pt.col=as.numeric(lcush$Type)+1,
       pt.pch=as.character(lcush$Type))

# miscellaneous other examples

tita <- get.tita()

mod.glm.sex <- glm(sex~., data=tita, family=binomial)
plotmo1(mod.glm.sex, pt.col=as.numeric(tita$pclass)+1, do.par=2)
empty.plot()

# glm binomial model: compare error band on plotmo versus manual plot

wasp <- read.table("darlingtonia.tab", header=T)
iorder = order(wasp$leafHeight, wasp$visited)
wasp = wasp[iorder,]
iseq = seq(1, nrow(wasp), by = 2)
dat = wasp[iseq, ]
newdat = wasp[iseq+1, ]

print(summary(dat))
set.seed(2026)
mod <- glm(visited ~ leafHeight, data=dat, family=binomial())
print(summary(mod))
set.seed(2026)
level95 = .95
plotmo(mod, do.par=FALSE, trace=2,
    level=level95, level.shade="lightgray",
    # nrug=TRUE,
    grid.col="lightgray",
    pt.col=ifelse(dat$visited, "pink", "lightgreen"), jitter=0.1,
    main="plotmo: visited vs leafHeight")

# wrong error bars (used prior to plotmo 3.7.0)
pred = predict(mod, newdata=newdat, type='response', se.fit=TRUE)
yhat = pred$fit
wrong_se = pred$se.fit
wrong_lwr = yhat - (2 * wrong_se)
wrong_upr = yhat + (2 * wrong_se)
points(newdat$leafHeight, wrong_lwr, pch=".", col="red")
points(newdat$leafHeight, wrong_upr, pch=".",, col="red")

# correct error bars (fixed in plotmo 3.7.0)
ilink = family(mod)$linkinv
pred_link = predict(mod, newdata=newdat, type="link", se.fit=TRUE)
q = qt((1 - level95) / 2, df = df.residual(mod), lower.tail = FALSE)
printf("darlingtonia: level95 %.4f q %.4f\n", level95, q)
lwr = ilink(pred_link$fit - (q * pred_link$se.fit))
upr = ilink(pred_link$fit + (q * pred_link$se.fit))

points(newdat$leafHeight, lwr, pch="o")
points(newdat$leafHeight, upr, pch="o")

legend(x="topleft",
       legend=c("correct (plotmo 3.7.0)", "incorrect (prior to plotmo 3.7.0)"),
       col=c(1,2),
       pch=sprintf("%s%s", "o", "."))

plot(newdat$leafHeight, yhat, type="l", main="manual: visited vs leafHeight",
     xlim=c(15,80), ylim=c(-.05,1.1))
points(newdat$leafHeight, lwr, pch="o")
points(newdat$leafHeight, upr, pch="o")

# tita[,4] is age, tita[,1] is pclass
printf("library(lars)\n")
library(lars)
set.seed(2015)
xmat <- as.matrix(tita[,c(2,5,6)])
mod.lars.xmat <- lars(xmat, tita[,4])
par(mfrow=c(2,2))
plot(mod.lars.xmat)
plotmo1(mod.lars.xmat, nresponse=4, do.par=F)
plotres(mod.lars.xmat, trace=0, nresponse=4)

if(0) { # TODO fails with R-3.4.2: object '.QP_qpgen2' not found
    printf("library(cosso)\n")
    library(cosso)
    set.seed(2016)
    cosso <- cosso(xmat,tita[,4],family="Gaussian")
    # TODO tell maintainer of cosso that you have to do this
    class(cosso) <- "cosso"
    set.seed(2016)
    plotmo1(cosso)
    set.seed(2016)
    plotres(cosso)
}
# examples from James, Witten, et al. ISLR book
# I tested all models in their scripts manually.
# All worked except for exceptions below.

printf("library(pls)\n")
library(pls)
printf("library(ISLR)\n")
library(ISLR)
Hitters=na.omit(Hitters)

set.seed(1)
x <- model.matrix(Salary~.,Hitters)[,-1]
y <- Hitters$Salary
train=sample(1:nrow(x), nrow(x)/2)
pcr.fit1=pcr(Salary~., data=Hitters,subset=train,scale=TRUE, validation="CV")
plotmo1(pcr.fit1, nresponse=10)

# set.seed(1)
# x <- model.matrix(Salary~.,Hitters)[,-1]
# y <- Hitters$Salary
# train=sample(1:nrow(x), nrow(x)/2)
# pcr.fit2=pcr(y~x,scale=TRUE,ncomp=7)
# # TODO following gives Error: predictions returned the wrong length (got 263 but expected 50)
# plotmo1(pcr.fit2, nresponse=5)

library(splines)
fit.lm2=lm(wage~bs(age,knots=c(25,40,60)),data=Wage)
par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),oma=c(0,0,4,0))
agelims=range(Wage$age)
age.grid=seq(from=agelims[1],to=agelims[2])
pred=predict(fit.lm2,newdata=list(age=age.grid),se=T)
plot(Wage$age,Wage$wage,col="gray", ylim=c(0,320))
lines(age.grid,pred$fit,lwd=2)
lines(age.grid,pred$fit+2*pred$se,lty="dashed")
lines(age.grid,pred$fit-2*pred$se,lty="dashed")
fit.lm2=lm(wage~bs(age,knots=c(25,40,60)),data=Wage,model=F) # TODO delete
plotmo1(fit.lm2, col.resp=2, do.par=F, level=.95, ylim=c(0,320),
        nrug=TRUE, caption="fit.lm2", ylab="wage")

fit.glm2 <- glm(I(wage>250)~poly(age,4),data=Wage,family=binomial)
par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),oma=c(0,0,4,0))
agelims=range(Wage$age)
age.grid=seq(from=agelims[1],to=agelims[2])
# their plot
preds=predict(fit.glm2,newdata=list(age=age.grid),se=T)
pfit=exp(preds$fit)/(1+exp(preds$fit))
se.bands.logit = cbind(preds$fit+2*preds$se.fit, preds$fit-2*preds$se.fit)
se.bands = exp(se.bands.logit)/(1+exp(se.bands.logit))
preds=predict(fit.glm2,newdata=list(age=age.grid),type="response",se=T)
plot(Wage$age,I(Wage$wage>250),xlim=agelims,type="n",ylim=c(0,.2))
points(jitter(Wage$age), I((Wage$wage>250)/5),cex=.5,pch="|",col="darkgrey")
lines(age.grid,pfit,lwd=2, col="blue")
matlines(age.grid,se.bands,lwd=1,col="blue",lty=3)
# plotmo plot, side by side
# TODO Warning: the level argument may not be properly supported on glm objects built with weights
plotmo1(fit.glm2, level=.95, degree1.col="blue", ylim=c(0,.2), do.par=FALSE, nrug=-1, caption="fit.glm2", ylab="I(wage > 250)")

# Test deparsing of the formula in plotmo.pairs.default
# TODO Height is included in the plots even though formula says -Height
Height2 <- trees$Height^2
a <- lm(Volume~(Girth*Height2)-Height, data=trees, x=TRUE, model=FALSE)
plotmo(a)

# test "the variable on the right side of the formula is a matrix or data.frame"
# TODO would like to solve this problem

options(warn=2)
data(gasoline, package="pls")
earth.octane <- earth(octane ~ NIR, data=gasoline)
print(summary(earth.octane)) # ok
plotres(earth.octane) # ok
expect.err(try(plotmo(earth.octane)), "the variable on the right side of the formula is a matrix or data.frame")
options(warn=1)

# TODO May 2020 'ElemStatLearn' is not available (for R version 4.0.0)
# library(ElemStatLearn)
# x <- mixture.example$x
# g <- mixture.example$y
# lm.mixture.example <- lm(g ~ x)
# options(warn=2)
# expect.err(try(plotmo(lm.mixture.example)), "the variable on the right side of the formula is a matrix or data.frame")
# options(warn=1)

# test variable names with $ are not supported

a <- earth(O3~ozone1$doy, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- earth(O3~ozone1$doy + temp, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- lm(O3~ozone1$doy, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

a <- lm(O3~ozone1$doy + temp, data=ozone1)
expect.err(try(plotmo(a)), "cannot get the original model predictors")

#--- test interaction of w1. and non w1 args -------------------------------------

par(mfrow=c(4,3), mar=c(3, 3, 4, 1), mgp=c(2, 0.6, 0))

mod78 <- earth(Volume ~ ., data = trees)
par(mfrow=c(3,4), mar=c(3, 3, 3, 1), mgp=c(2, 0.6, 0))

# multiple which, earth model
plotres(mod78, cex.main=1,
        ylim=c(-.5, .8),    xlim=c(-2, 7),    col=2:3, do.par=FALSE,
        w1.main=c("ylim=c(-.5, .8)\nxlim=c(-2, 7) col=2:3"))

# multiple which, earth model
plotres(mod78, cex.main=.7,
        w1.ylim=c(-.5, .8), w1.xlim=c(-2, 7), col=2:3, do.par=FALSE,
        ylim=c(-10,10), xlim=c(-30, 100),
        w1.main=c("w1.ylim=c(-.5, .8) w1.xlim=c(-2, 7)\nylim=c(-10,10), xlim=c(-30, 100)"))
par(org.par)

par(mfrow=c(3,4), mar=c(3, 3, 3, 1), mgp=c(2, 0.6, 0))

# which=1, earth model

plotres(mod78, which=1, cex.main=.8,
        col=2:3,
        main="which=1, no other ylim args",
        w1.main="which=1, no other ylim args")

plotres(mod78, which=1, cex.main=.8,
        col=2:3, w1.ylim=c(.3,.98), w1.xlim=c(-2, 7),
        main="w1.ylim=c(.3,.98)\nw1.xlim=c(-2, 7)")

plotres(mod78, which=1, cex.main=.8,
        col=2:3, ylim=c(.3,.98),    xlim=c(-2, 7),
        main="ylim=c(.3,.98)\nxlim=c(-2, 7)")  # ylim gets passed to modsel

plotres(mod78, which=1, cex.main=.75,
        col=2:3, w1.ylim=c(.3,.98), ylim=c(-.5,.5),
        w1.xlim=c(-2, 7), xlim=c(-90, 90),
        main="w1.ylim=c(.3,.98), ylim=c(-.5,.5)\nw1.xlim=c(-2, 7), xlim=c(-90, 90)") # ignore ylim

# which=3, earth model
plotres(mod78, which=3, cex.main=1,
        col=2:3,
        main="which=3, no other ylim args")

plotres(mod78, which=3, cex.main=1,
        col=2:3, w1.ylim=c(.3,.98), w1.xlim=c(-2, 7),
        main="w1.ylim=c(.3,.98)\nw1.xlim=c(-2, 7)") # not usual, ignore w1.ylim

plotres(mod78, which=3, cex.main=1,
        col=2:3, ylim=c(-10,10), xlim=c(-90,90),
         main="which=3, ylim=c(-10,10)\nxlim=c(-90,90)")

plotres(mod78, which=3, cex.main=1,
        col=2:3, w1.ylim=c(.3,.98), ylim=c(-10,10), w1.xlim=c(-2, 7), xlim=c(-90,90),
        main="w1.ylim=c(.3,.98) ylim=c(-10,10)\nw1.xlim=c(-2, 7), xlim=c(-90,90)")

par(org.par)

nullarg <- NULL
expect.err(try(plotmo(nullarg)),   "argument 'nullarg' is NULL")
expect.err(try(plotmo(NULL)),      "argument 'NULL' is NULL")
expect.err(try(plotmo(0)),         "'0' is not an S3 model")
expect.err(try(plotmo(list(1,2))), "'list(1, 2)' is a plain list, not an S3 model")
expect.err(try(plotmo(list(1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0))),
               "object is a plain list, not an S3 model")

source("test.epilog.R")
