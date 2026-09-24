# test.mlr.R: test the "mlr" package with plotmo and plotres
#
# TODO mlr is in maintenance mode, add mlr3 support to plotmo?
# TODO generally, plotres residuals for WrappedModel prob models aren't right

source("test.prolog.R")
library(mlr)
library(plotmo)
library(rpart.plot)
library(earth)
# TODO following function is temporary until mlr package is updated
train.with.call <- function(learner, task, subset=NULL, weights=NULL)
{
    retval <- train(learner, task, subset, weights)
    retval$call <- match.call()
    retval
}

cat("==simple one variable regression model with earth ===============================\n")

data(trees)
trees1 <- trees[,c("Volume", "Girth")]

task <- makeRegrTask(data=trees1, target="Volume")
lrn <- makeLearner("regr.earth", degree=2)
regr.earth.with.call = train.with.call(lrn, task)
regr.earth = train(lrn, task)
earth <- earth(Volume~., data=trees1, degree=2)

# SHOWCALL is just a testing thing, so we can see who created the plot on the plot itself
plotres(regr.earth.with.call, SHOWCALL=TRUE)
plotres(regr.earth$learner.model, SHOWCALL=TRUE)
plotres(earth, SHOWCALL=TRUE)

plotmo(regr.earth.with.call, trace=1, SHOWCALL=TRUE)
plotmo(regr.earth$learner.model, trace=1, SHOWCALL=TRUE)
plotmo(earth, trace=1, SHOWCALL=TRUE)

# compare partial dependence plots from mlr and plotmo packages
set.seed(2018)
plotmo(earth, pmethod="partdep", SHOWCALL=TRUE, col=2, pt.col="darkgray", grid.col="lightgray")
set.seed(2018)
pd <- generatePartialDependenceData(regr.earth, task, "Girth", n=c(50, NA))
print(plotPartialDependence(pd, data = getTaskData(task)))

cat("==test error handling if original data is messed up===========================\n")

par(mfrow=c(4,2), mar=c(1.5,2.5,4,1), oma=c(0,0,0,0))
colnames(trees1) <- c("nonesuch", "Volume")
plotmo(regr.earth$learner.model, do.par=0, degree1=1, degree2=0, main='colnames(trees1) <- c("nonesuch", "Volume")')
plotmo(regr.earth.with.call, do.par=0, degree1=1, degree2=0)
par(org.par)
expect.err(try(plotmo(earth, degree1=1, degree2=0)), "cannot get the original model predictors")

cat("==regression model with randomForest (binary response)============================\n")

library(randomForest)
library(earth) # for etitanic data
data(etitanic)
set.seed(2018)
# use a logical subset (since we test for numeric subset elsewhere)
# use a small subset so we can see easily if subset is applied or ignored in plots
train.subset <- rnorm(nrow(etitanic)) > 1 # 166 cases ((16% of 1046 cases))
printf("sum(train.subset) %g (%.0f%% of %g cases)\n", sum(train.subset),
    100 * sum(train.subset) / nrow(etitanic), nrow(etitanic))
task.regr.rf <- makeRegrTask(data=etitanic, target="survived")
lrn.regr.rf = makeLearner("regr.randomForest")
set.seed(2018)
regr.rf.with.call = train.with.call(lrn.regr.rf, task.regr.rf, subset=train.subset)
set.seed(2018)
rf <- randomForest(survived~., data=etitanic, subset=train.subset)
# sanity check that the models are identical
stopifnot(identical(predict(regr.rf.with.call$learner.model), predict(rf)))

plotres(regr.rf.with.call, info=TRUE, SHOWCALL=TRUE)
# plotres(regr.rf$learner.model, info=TRUE, SHOWCALL=TRUE) # Error: no formula in getCall(object)
plotres(rf, info=TRUE, SHOWCALL=TRUE)

set.seed(2018) # for repeatable jitter in points (specified with pt.col)
plotmo(regr.rf.with.call, pt.col=2, SHOWCALL=TRUE)
# plotmo(regr.rf$learner.model, trace=1, SHOWCALL=TRUE) # Error: no formula in getCall(object)
set.seed(2018)
plotmo(rf, pt.col=2, SHOWCALL=TRUE)

# compare partial dependence plots
set.seed(2018)
plotmo(regr.rf.with.call, degree1="age", degree2=0, pmethod="partdep",
       grid.col="gray", col=2, pt.col="darkgray", SHOWCALL=TRUE)
# function from randomForest package
set.seed(2018)
partialPlot(rf, pred.data=etitanic[train.subset,], x.var="age", n.pt=50, ylim=c(0, 1))
grid()
# function from mlr package
set.seed(2018)
pd <- generatePartialDependenceData(regr.rf.with.call, task.regr.rf, "age", n=c(50, NA))
print(plotPartialDependence(pd, data = getTaskData(task.regr.rf)))

plotmo(regr.rf.with.call, degree1="pclass", degree2=0, pmethod="partdep", SHOWCALL=TRUE)
set.seed(2018)
# function from randomForest package
set.seed(2018)
partialPlot(rf, pred.data=etitanic[train.subset,], x.var="pclass", n.pt=50, ylim=c(0, 1))
grid()
# TODO following fails
pd <- generatePartialDependenceData(regr.rf.with.call, task.regr.rf, "pclass", n=c(50, NA))
try(print(plotPartialDependence(pd, data = getTaskData(task.regr.rf)))) # Error: Discrete value supplied to continuous scale

cat("==classification model with randomForest (binary response)======================\n")

set.seed(2018)
library(earth) # for etitanic data
data(etitanic)
etit <- etitanic
etit$survived <- factor(etit$survived, labels=c("notsurvived", "survived"))

task.classif.rf <- makeClassifTask(data=etit, target="survived")
lrn.classif.rf <- makeLearner("classif.randomForest", predict.type="prob")
set.seed(2018)
classif.rf.with.call <- train.with.call(lrn.classif.rf, task.classif.rf, , subset=train.subset)
set.seed(2018)
rf <- randomForest(survived~., data=etit, method="class", subset=train.subset)
# sanity check that the models are identical
stopifnot(identical(predict(classif.rf.with.call$learner.model), predict(rf)))

# TODO following causes Error: classif.earth: Setting parameter glm without available description object
# lrn <- makeLearner("classif.earth", degree=2, glm=list(family=binomial))

# TODO residuals on WrappedModel don't match direct call to rf model
set.seed(2018) # for repeatable jitter
plotres(classif.rf.with.call, nresponse="prob.survived", SHOWCALL=TRUE, jitter=2)
set.seed(2018)
plotres(classif.rf.with.call$learner.model, type="prob", SHOWCALL=TRUE, jitter=2)
set.seed(2018)
plotres(rf, type="prob", SHOWCALL=TRUE, jitter=2)

options(warn=2) # treat warnings as errors
expect.err(try(plotmo(classif.rf.with.call)), "Defaulting to nresponse=1, see above messages")
options(warn=1)
set.seed(2018) # for repeatable jitter
plotmo(classif.rf.with.call,               SHOWCALL=TRUE, nresponse="prob.survived", pt.col=2, trace=2)
set.seed(2018)
plotmo(classif.rf.with.call$learner.model, SHOWCALL=TRUE, type="prob", pt.col=2)
set.seed(2018)
# note that in the following, get.y.shift.scale (in plotmo code) rescales the plotted y to 0..1
plotmo(rf,                                 SHOWCALL=TRUE, type="prob", pt.col="gray")
set.seed(2018)
# in following graph, note that get.y.shift.scale doesn't rescale the plotted y because ylim=c(0,2)
plotmo(rf,                                 SHOWCALL=TRUE, type="prob", ylim=c(0,2), pt.col="gray")

# compare partial dependence plots
set.seed(2018)
plotmo(rf, type="prob", degree1="pclass", degree2=0, pmethod="partdep", pt.col=2, SHOWCALL=TRUE)
set.seed(2018)
plotmo(rf,              degree1="pclass", degree2=0, pmethod="partdep", pt.col=2, SHOWCALL=TRUE)
set.seed(2018)
# TODO following fails
pd <- generatePartialDependenceData(classif.rf.with.call, task.classif.rf, "pclass", n=c(50, NA))
try(print(plotPartialDependence(pd, data = getTaskData(task.classif.rf)))) # Error: Discrete value supplied to continuous scale

plotmo(rf, type="prob", nresponse="notsurvived", degree1="age", degree2=0,
       pmethod="partdep", ylim=c(.3,.75), nrug=TRUE, grid.col="gray") # looks plausible
set.seed(2018)
pd <- generatePartialDependenceData(classif.rf.with.call, task.classif.rf, "age", n=c(50, NA))
print(plotPartialDependence(pd, data = getTaskData(task.classif.rf)))

cat("==examples from plotmo-notes.pdf ===============================================\n")

#-- Regression model with mlr -------------------------------------------

library(mlr)
library(plotmo)
lrn <- makeLearner("regr.svm")
fit1.with.call <- train.with.call(lrn, bh.task)
fit1 <- train(lrn, bh.task)

# generate partial dependence plots for all variables
# we use "apartdep" and not "partdep" to save testing time
plotmo(fit1.with.call, pmethod="apartdep")
plotmo(fit1$learner.model, pmethod="apartdep")

# generate partial dependence plot for just "lstat"
set.seed(2018) # so slight jitter on pt.col points in plotmo doesn't change across test runs
plotmo(fit1.with.call,
       degree1="lstat",           # what predictor to plot
       degree2=0,                 # no interaction plots
       pmethod="partdep",         # generate partial dependence plot
       pt.col=2, grid.col="gray", # optional bells and whistles
       nrug=TRUE)                 # rug ticks along the bottom
set.seed(2018) # so slight jitter on pt.col points in plotmo doesn't change across test runs
plotmo(fit1$learner.model,
       degree1="lstat",           # what predictor to plot
       degree2=0,                 # no interaction plots
       pmethod="partdep",         # generate partial dependence plot
       pt.col=2, grid.col="gray", # optional bells and whistles
       nrug=TRUE)                 # rug ticks along the bottom

# compare to the function provided by the mlr package
set.seed(2018)
pd <- generatePartialDependenceData(fit1, bh.task, "lstat", n=c(50, NA))
print(plotPartialDependence(pd, data = getTaskData(bh.task)))
# # TODO following fails: Error: Discrete value supplied to continuous scale
# pd <- generatePartialDependenceData(fit1, bh.task, "chas", n=c(50, NA))
# plotPartialDependence(pd, data = getTaskData(bh.task))

#-- Classification model with mlr ---------------------------------------

lrn.classif.rpart <- makeLearner("classif.rpart", predict.type = "prob", minsplit = 10)
fit2.with.call <- train.with.call(lrn.classif.rpart, iris.task)
fit2 <- train(lrn.classif.rpart, iris.task)

# generate partial dependence plots for all variables
# TODO plotmo can plot the response for only one class at a time
plotmo(fit2.with.call,
       nresponse="prob.virginica", # what response to plot
       # type="prob",              # type gets passed to predict.rpart
       pmethod="apartdep")         # generate partial dependence plot

plotmo(fit2$learner.model,
       nresponse="virginica",     # what response to plot
       type="prob",               # type gets passed to predict.rpart
       pmethod="apartdep")        # generate partial dependence plot

# generate partial dependence plot for just "Petal.Length"
plotmo(fit2.with.call,
       degree1="Petal.Length",    # what predictor to plot
       degree2=0,                 # no interaction plots
       nresponse="prob.virginica",     # what response to plot
       # type="prob",               # type gets passed to predict.rpart
       pmethod="apartdep")        # generate partial dependence plot

plotmo(fit2$learner.model,
       degree1="Petal.Length",    # what predictor to plot
       degree2=0,                 # no interaction plots
       nresponse="virginica",     # what response to plot
       type="prob",               # type gets passed to predict.rpart
       pmethod="apartdep")        # generate partial dependence plot

# compare to the function provided by the mlr package
set.seed(2018)
pd <- generatePartialDependenceData(fit2, iris.task, "Petal.Length", n=c(50, NA))
print(plotPartialDependence(pd, data = getTaskData(iris.task)))

cat("==lda example from mlr documentation, and plotmo error handling =================\n")

set.seed(2018)
data(iris)
task.lda <- makeClassifTask(data=iris, target="Species")
lrn.lda <- makeLearner("classif.lda")
n <- nrow(iris)
train.set <- sample(n, size=2/3*n)
test.set <- setdiff(1:n, train.set)
classif.lda.with.call <- train.with.call(lrn.lda, task.lda, subset=train.set)
classif.lda <- train(lrn.lda, task.lda, subset=train.set)
iris1 <- iris[train.set, ]
library(MASS)
lda <- lda(Species~., data=iris1)

# expect.err(try(plotres(classif.lda.with.call)), "plotres does not (yet) support type=\"class\" for \"lda\" objects")
expect.err(try(plotres(classif.lda$learner.model)), "plotres does not (yet) support type=\"class\" for \"lda\" objects")

options(warn=2) # treat warnings as errors
# expect.err(try(plotres(classif.lda.with.call, type="response")), "predict.lda returned multiple columns (see above) but nresponse is not specified")
expect.err(try(plotres(classif.lda$learner.model, type="response")), "Defaulting to nresponse=1, see above messages")
options(warn=1)

expect.err(try(plotres(classif.lda.with.call, type="response", nresponse="nonesuch")), "nresponse=\"nonesuch\" is not allowed")
expect.err(try(plotres(classif.lda$learner.model, type="response", nresponse="nonesuch")), "nresponse=\"nonesuch\" is not allowed")

expect.err(try(plotres(classif.lda.with.call, type="response", nresponse=0)), "nresponse=0 but it should be at least 1")
expect.err(try(plotres(classif.lda$learner.model, type="response", nresponse=0)), "nresponse=0 but it should be at least 1")

expect.err(try(plotres(classif.lda.with.call, type="response", nresponse=99)), "nresponse is 99 but the number of columns is only 1")
expect.err(try(plotres(classif.lda$learner.model, type="response", nresponse=99)), "nresponse is 99 but the number of columns is only 2")

expect.err(try(plotmo(classif.lda)), "getCall(classif.lda) failed")

expect.err(try(plotres(classif.lda)), "getCall(classif.lda) failed")

# TODO residuals don't match
plotres(classif.lda.with.call,     SHOWCALL=TRUE, type="response")
plotres(classif.lda$learner.model, SHOWCALL=TRUE, type="response", nresponse="LD2")
plotres(lda,                       SHOWCALL=TRUE, type="response", nresponse="LD2")

plotmo(classif.lda.with.call,     SHOWCALL=TRUE)
plotmo(classif.lda$learner.model, SHOWCALL=TRUE)
plotmo(lda,                       SHOWCALL=TRUE)

# # TODO plotPartialDependence and plotmo graphs below don't match
# pd <- generatePartialDependenceData(classif.lda, task.lda, "Petal.Width", n=c(50, NA)) # TODO generates warnings
# print(plotPartialDependence(pd, data = getTaskData(task.lda)))
plotmo(classif.lda.with.call, degree1="Petal.Width", degree2=0, pmethod="partdep", do.par=FALSE)

plotmo(classif.lda.with.call,     SHOWCALL=TRUE, all2=TRUE, type="response")
plotmo(classif.lda$learner.model, SHOWCALL=TRUE, all2=TRUE, type="class")
plotmo(lda,                       SHOWCALL=TRUE, all2=TRUE, type="class")

plotmo(classif.lda$learner.model, SHOWCALL=TRUE, all2=TRUE, type="response", nresponse="LD1")
plotmo(lda,                       SHOWCALL=TRUE, all2=TRUE, type="response", nresponse="LD1")

cat("==test recursive call to plotmo_prolog for learner.model===============\n")

set.seed(2018)
n <- 100
data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = rnorm(n),
    x5 = rnorm(n),
    x6 = rnorm(n),
    x7 = rnorm(n),
    x8 = rnorm(n),
    x9 = rnorm(n))

data$y <- sin(data$x3) + sin(data$x4) + 2 * cos(data$x5)

set.seed(2018)
library(gbm)
# reference model
gbm = gbm(y~., data=data, n.trees=300)
plotmo(gbm, trace=-1, SHOWCALL=TRUE)

set.seed(2018)
task <- makeRegrTask(data=data, target="y")
lrn <- makeLearner("regr.gbm", n.trees=300, keep.data=TRUE)
regr.gbm = train.with.call(lrn, task)
plotmo(regr.gbm, trace=-1, SHOWCALL=TRUE)

set.seed(2018)
lrn <- makeLearner("regr.gbm", n.trees=300)
regr.gbm.nokeepdata = train.with.call(lrn, task)
# expect message: use keep.data=TRUE in the call to gbm (cannot determine the variable importances)
plotmo(regr.gbm.nokeepdata, trace=1, SHOWCALL=TRUE)

plotres(regr.gbm, SHOWCALL=TRUE)

cat("==example from makeClassificationViaRegressionWrapper help page ===============\n")
# this tests that plotmo.prolog can access the learner.model at object$learner.model$next.model$learner.model

set.seed(2018)
lrn = makeLearner("regr.rpart")
lrn = makeClassificationViaRegressionWrapper(lrn)
ClassificationViaRegression = train.with.call(lrn, sonar.task, subset = 1:140)
plotmo(ClassificationViaRegression, SHOWCALL=TRUE)

source("test.epilog.R")
