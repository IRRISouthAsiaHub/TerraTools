# test.gbm.R: gbm tests for plotmo and plotres

source("test.prolog.R")
library(gbm)
library(rpart.plot) # for ptitanic, want data with NAs for testing
library(plotmo)
data(ptitanic)

cat("--- distribution=\"gaussian\", formula interface ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ] # small data for fast test
set.seed(2016)
# # TODO bug in gbm: following causes error: survived is not of type numeric, ordered, or factor
# ptit$survived <- ptit$survived == "survived"
ptit <- ptit[!is.na(ptit$age), ]
train.frac <- .8
set.seed(2016)
gbm.gaussian <- gbm(age~., data=ptit, train.frac=train.frac,
                   distribution="gaussian",
                   n.trees=50, shrinkage=.1, keep.data=FALSE)
expect.err(try(plotres(gbm.gaussian)), "use keep.data=TRUE in the call to gbm")
set.seed(2016)
gbm.gaussian <- gbm(age~., data=ptit, train.frac=train.frac,
                   distribution="gaussian",
                   n.trees=50, shrinkage=.1)
par(mfrow=c(2,2), mar=c(3,3,4,1))
w1 <- plotres(gbm.gaussian, which=1, do.par=FALSE, w1.smooth=TRUE,
        w1.main="gbm.gaussian")
cat("w1 plot for gbm.gaussian returned (w1.smooth=TRUE):\n")
print(w1)
plot(0, 0) # dummy plot
w3 <- plotres(gbm.gaussian, which=3, do.par=FALSE, info=TRUE,
        smooth.col=0, col=ptit$sex, # ylim=c(-40,40),
        wmain="nresponse=1")
# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$age[iused]
n.trees <- plotmo:::gbm.n.trees(gbm.gaussian)
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm.gaussian, type="response", n.trees=n.trees)
yhat <- predict(gbm.gaussian, newdata=ptit, type="response", n.trees=n.trees)
yhat <- yhat[iused]
plot(yhat, y - yhat,
     col=ptit$sex[iused], main="manual gaussian residuals",
     pch=20, ylim=c(-40,40))
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
par(org.par)

w1 <- plotres(gbm.gaussian, predict.n.trees=13, w1.grid.col=1, trace=1, SHOWCALL=TRUE,
              w1.smooth=TRUE,
              w1.main="predict.n.trees=13 w1.grid.col=1")
cat("second w1 plot for gbm.gaussian returned (w1.smooth=TRUE):\n")
print(w1)
plotmo(gbm.gaussian, trace=-1, SHOWCALL=TRUE)
# plotmo(gbm.gaussian, trace=-1, all1=TRUE, SHOWCALL=TRUE)
# plotmo(gbm.gaussian, trace=-1, all2=TRUE, SHOWCALL=TRUE)

# test color argument
par(mfrow=c(2,2), mar=c(3,3,4,1))
plotres(gbm.gaussian, which=1)
title("test color argument")
plotres(gbm.gaussian, which=1, w1.col=c(1,2,3,0))
plotres(gbm.gaussian, which=1, w1.col=c(1,0,0,4), w1.legend.x=40, w1.legend.y=.3)
plotres(gbm.gaussian, which=1, w1.col=c(2,3,4,1), w1.legend.x="topright")
par(org.par)

par(mfrow=c(2,2), mar=c(3,3,4,1))
plot_gbm(gbm.gaussian)
title("test plot_gbm")
w1 <- plot_gbm(gbm.gaussian, col=c(1,2,3,0), grid.col=1, smooth=TRUE,
               main="col=c(1,2,3,0), grid.col=1")
cat("third w1 plot for gbm.gaussian returned (smooth=TRUE):\n")
print(w1)
par(org.par)

# test xlim and ylim
par(mfrow=c(2,3), mar=c(3,3,4,1))
plot_gbm(gbm.gaussian,                  main="test xlim and ylim default")
plot_gbm(gbm.gaussian, ylim=NULL,       main="ylim=NULL")
plot_gbm(gbm.gaussian, xlim=c(5, 50),    main="xlim=c(5, 50)")
plot_gbm(gbm.gaussian, ylim=c(100, 250), main="ylim=c(100, 250)")
plot_gbm(gbm.gaussian, xlim=c(10, 25),
                       ylim=c(150, 170),  main="xlim=c(10, 25), ylim=c(150, 170)")
plot_gbm(gbm.gaussian, xlim=c(-10, 40), ylim=c(-10, 300), legend.x=NA,
                                        main="xlim=c(-10, 40), ylim=c(-10, 300)\nlegend.x=NA")
par(org.par)

# test the smooth argument
par(mfrow=c(3,3), mar=c(3,3,4,1))
imin <- plot_gbm(gbm.gaussian,                    main="smooth=default")
imin.default <- imin
cat("smooth=default imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(1,0,0,0), main="smooth=c(1,0,0,0)")
cat("smooth=c(1,0,0,0) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(0,1,0,0), main="smooth=c(0,1,0,0)")
cat("smooth=c(0,1,0,0) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(0,0,1,0), main="smooth=c(0,0,1,0)")
cat("smooth=c(0,0,1,0) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(0,0,0,1), main="smooth=c(0,0,0,1)\nsame as default")
cat("smooth=c(0,0,0,1) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(0,0,0,0), main="smooth=c(0,0,0,0)")
cat("smooth=c(0,0,0,0) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=c(0,0,1,1), main="smooth=c(0,0,1,1)")
cat("smooth=c(0,0,1,1) imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")

imin <- plot_gbm(gbm.gaussian, smooth=1,          main="smooth=1") # gets recycled
cat("smooth=1          imin=c(", imin[1], ",", imin[2], ",", imin[3], ",", imin[4], ")\n", sep="")
imin.smooth <- imin

imin.noplot <- plot_gbm(gbm.gaussian, col=0) # will not be plotted
print(imin.default)
print(imin.noplot)
stopifnot(identical(imin.default, imin.noplot))

imin.noplot <- plot_gbm(gbm.gaussian, col=0, smooth=1) # will not be plotted
print(imin.smooth)
print(imin.noplot)
stopifnot(identical(imin.smooth, imin.noplot))

par(org.par)

cat("--- distribution=\"gaussian\", glm.fit interface ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ]
set.seed(2016)
ptit <- ptit[!is.na(ptit$age), ]
train.frac <- .8
set.seed(2016)
gbm.gaussian.fit <- gbm.fit(ptit[,-4], ptit[,4], nTrain=floor(train.frac * nrow(ptit)),
                   distribution="gaussian", verbose=FALSE,
                   n.trees=50, shrinkage=.1)
par(mfrow=c(2,2), mar=c(3,3,4,1))
w1 <- plotres(gbm.gaussian.fit, which=1, do.par=FALSE, w1.smooth=TRUE,
        w1.main="gbm.gaussian.fit")

cat("w1 plot for gbm.gaussian.fit returned (w1.smooth=TRUE):\n")
print(w1)

plot(0, 0) # dummy plot

w3 <- plotres(gbm.gaussian.fit, which=3, do.par=FALSE, info=TRUE, trace=0,
        smooth.col=0, col=ptit$sex, # ylim=c(-40,40),
        wmain="nresponse=1")

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y.fit <- ptit$age[iused]
n.trees <- plotmo:::gbm.n.trees(gbm.gaussian.fit)
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat.fit <- predict(gbm.gaussian.fit, type="response", n.trees=n.trees)
yhat.fit <- predict(gbm.gaussian.fit, newdata=ptit[,-4], type="response", n.trees=n.trees)
yhat.fit <- yhat.fit[iused]
# plot(yhat.fit, y.fit - yhat.fit,
#      col=ptit$sex[iused], main="manual gaussian residuals\n(TODO gbm.fit don't match)",
#      pch=20, ylim=c(-40,40))
# abline(h=0, col="gray")
# --- TODO known issue, these fail ---
# compare to formual interface
# stopifnot(all(yhat.fit == yhat))
stopifnot(all(y.fit == y))
# # sanity check
# stopifnot(all(yhat.fit == w3$x))
# stopifnot(all(y.fit - yhat.fit == w3$y.fit))
plotmo(gbm.gaussian.fit, trace=-1, SHOWCALL=TRUE)
par(org.par)

cat("--- distribution=\"laplace\" ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ]
ptit <- ptit[!is.na(ptit$age), ]
ptit$survived <- ptit$parch <- ptit$sex <- NULL
train.frac <- .8
set.seed(2016)
gbm.laplace <- gbm(age~., data=ptit, train.frac=train.frac,
                   distribution="laplace",
                   n.trees=100, shrinkage=.1)
par(mfrow=c(2,2), mar=c(3,3,4,1))
w1 <- plotres(gbm.laplace, which=1:2, do.par=FALSE, w1.smooth=TRUE,
        w1.main="gbm.laplace")

cat("w1 plot for gbm.laplace returned (w1.smooth=TRUE):\n")
print(w1)

w3 <- plotres(gbm.laplace, which=3, do.par=FALSE, info=TRUE)

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$age[iused]
n.trees <- plotmo:::gbm.n.trees(gbm.laplace)
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm.laplace, type="response", n.trees=n.trees)
yhat <- predict(gbm.laplace, newdata=ptit, type="response", n.trees=n.trees)
yhat <- yhat[iused]
plot(yhat, y - yhat,
     main="manual laplace residuals",
     pch=20, ylim=c(-40,40))
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
plotmo(gbm.laplace, trace=-1, SHOWCALL=TRUE)
par(org.par)

# # TODO commented out because gives random slightly different results per invocation
# cat("--- distribution=\"tdist\" ----------------------------------\n")
#
# set.seed(2016)
# ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ]
# ptit <- ptit[!is.na(ptit$age), ]
# ptit$survived <- ptit$parch <- ptit$sex <- NULL
# train.frac <- .8
# set.seed(2016)
# gbm.tdist <- gbm(age~., data=ptit, train.frac=train.frac,
#                    distribution="tdist",
#                    n.trees=100, shrinkage=.1)
# par(mfrow=c(2,2), mar=c(3,3,4,1))
# set.seed(2016)
# w1 <- plotres(gbm.tdist, which=1:2, do.par=FALSE,
#         w1.main="gbm.tdist")
#
# cat("w1 plot for gbm.tdist returned (w1.smooth=default):\n")
# print(w1)
#
# set.seed(2016)
# w3 <- plotres(gbm.tdist, which=3, do.par=FALSE, info=TRUE)
#
# # compare to manual residuals
# iused <- 1:(train.frac * nrow(ptit))
# y <- ptit$age[iused]
# n.trees <- plotmo:::gbm.n.trees(gbm.tdist)
# # TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# # yhat <- predict(gbm.tdist, type="response", n.trees=n.trees)
# yhat <- predict(gbm.tdist, newdata=ptit, type="response", n.trees=n.trees)
# yhat <- yhat[iused]
# plot(yhat, y - yhat,
#      main="manual tdist residuals",
#      pch=20, ylim=c(-40,40))
# abline(h=0, col="gray")
# stopifnot(all(yhat == w3$x))
# stopifnot(all(y - yhat == w3$y))
# plotmo(gbm.tdist, trace=-1, SHOWCALL=TRUE)
# par(org.par)

cat("--- distribution=\"bernoulli\" ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=80), ]
ptit$survived <- as.numeric(ptit$survived == "survived")
temp <- ptit$pclass # put pclass at the end so can check ordering of importances
ptit$pclass <- NULL
ptit$pclass <- factor(as.numeric(temp), labels=c("first", "second", "third"))
train.frac <- .9
set.seed(2016)
gbm.bernoulli <- gbm(survived~., data=ptit, train.frac=train.frac,
                     distribution="bernoulli",
                     n.trees=100, shrinkage=.1, cv.folds=3)
par(mfrow=c(2,2))
par(mar=c(3.5, 3, 2, 0.5))  # small margins and text to pack figs in
par(mgp=c(1.5, .4, 0))      # squash axis annotations
w1 <- plotres(gbm.bernoulli, which=c(1,4),
        col=ptit$survived+2, trace=0, do.par=FALSE,
        w1.main="gbm.bernoulli")
cat("w1 plot for gbm.bernoulli with cv.folds=3 returned:\n")
print(w1)

w3 <- plotres(gbm.bernoulli, which=3, predict.n.trees=40,
        ylim=c(-.6, 1), xlim=c(.1, .6),
        col=ptit$sex, trace=0, do.par=FALSE, smooth.col=0)

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$survived[iused]
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm.bernoulli, type="response", n.trees=40)
yhat <- predict(gbm.bernoulli, newdata=ptit, type="response", n.trees=40)
yhat <- yhat[iused]
plot(yhat, y - yhat, col=ptit$sex,
     main="manual bernoulli residuals", pch=20, cex=1,
     ylim=c(-.6, 1), xlim=c(.1, .6))
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
par(org.par)

plotmo(gbm.bernoulli, do.par=2)
print(summary(gbm.bernoulli)) # will also plot
par(org.par)

cat("--- distribution=\"huberized\" ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=100), ]
ptit$survived <- as.numeric(ptit$survived == "survived")
ptit$sibsp <- ptit$parch <- ptit$pclass <- NULL
train.frac <- 1
set.seed(2016)
gbm.huberized <- gbm(survived~., data=ptit, train.frac=train.frac,
                     distribution="huberized",
                     n.trees=200, shrinkage=.1)
par(mfrow=c(2,2))
par(mar=c(3.5, 3, 2, 0.5))  # small margins and text to pack figs in
par(mgp=c(1.5, .4, 0))      # squash axis annotations
w1 <- plotres(gbm.huberized, which=c(1,4),
        col=ptit$survived+2, trace=0, do.par=FALSE,
        w1.main="gbm.huberized")
cat("w1 plot for gbm.huberized returned (smooth=default):\n")
print(w1)

# TODO huberized residuals look weird
w3 <- plotres(gbm.huberized, which=3, predict.n.trees=40,
        col=ptit$sex, trace=0, do.par=FALSE, smooth.col=0)

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$survived[iused]
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm.huberized, type="response", n.trees=40)
yhat <- predict(gbm.huberized, newdata=ptit, type="response", n.trees=40)
yhat <- yhat[iused]
plot(yhat, y - yhat, col=ptit$sex, ylim=c(-2.5, 2.5),
     main="manual huberized residuals", pch=20)
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
par(org.par)

plotmo(gbm.huberized, do.par=2)
print(summary(gbm.huberized)) # will also plot
par(org.par)

cat("--- distribution=\"adaboost\" ----------------------------------\n")

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=100), ]
ptit$survived <- as.numeric(ptit$survived == "survived")
ptit$sibsp <- ptit$parch <- ptit$pclass <- NULL
train.frac <- .8
set.seed(2016)
gbm.adaboost <- gbm(survived~., data=ptit, train.frac=train.frac,
                     distribution="adaboost",
                     n.trees=150, shrinkage=.01)
par(mfrow=c(2,2))
par(mar=c(3.5, 3, 2, 0.5))  # small margins and text to pack figs in
par(mgp=c(1.5, .4, 0))      # squash axis annotations
w1 <- plotres(gbm.adaboost, which=c(1,4),
        col=ptit$survived+2, trace=0, do.par=FALSE,
        w1.main="gbm.adaboost")
cat("w1 plot for gbm.adaboost returned (smooth=default):\n")
print(w1)

w3 <- plotres(gbm.adaboost, which=3, predict.n.trees=40,
        col=ptit$sex, trace=0, do.par=FALSE, smooth.col=0)

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$survived[iused]
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm.adaboost, type="response", n.trees=40)
yhat <- predict(gbm.adaboost, newdata=ptit, type="response", n.trees=40)
yhat <- yhat[iused]
plot(yhat, y - yhat, col=ptit$sex,
     main="manual adaboost residuals", pch=20)
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
par(org.par)

plotmo(gbm.adaboost, do.par=2)
print(summary(gbm.adaboost)) # will also plot
par(org.par)

# test gbm multinomial model, also test very small number of trees in plot_gbm

data(iris)
set.seed(2016)
gbm.iris <- gbm(Species~., data=iris, distribution="multinomial", n.tree=5)
expect.err(try(plotres(gbm.iris)),
           "gbm distribution=\"multinomial\" is not yet supported")
expect.err(try(plotmo(gbm.iris)),
           "gbm distribution=\"multinomial\" is not yet supported")
plot_gbm(gbm.iris)

# TODO following fails in the new version of gbm (version 2.2)
#   (distribution "multinomial" is no longer supported)
#
# cat("--- distribution=\"multinomial\" ----------------------------------\n")
#
# set.seed(2016)
# ptit <- ptitanic[sample(1:nrow(ptitanic), size=500), ]
# set.seed(2016)
# gbm.multinomial <- gbm(pclass~.,
#                        data=ptit, train.frac=.7,
#                        distribution="multinomial",
#                        n.trees=100, shrinkage=.1)
#
# w1 <- plot_gbm(gbm.multinomial, main="gbm.multinomial", smooth=T)
# cat("plot_gbm for gbm.multinomial returned (smooth=TRUE):\n")
# print(w1)
#
# expect.err(try(plotres(gbm.multinomial)),
#            "gbm distribution=\"multinomial\" is not yet supported")
#
# expect.err(try(plotmo(gbm.multinomial)),
#            "gbm distribution=\"multinomial\" is not yet supported")

# cat("--- gbmt distribution=\"Gaussian\", formula interface ----------------------------------\n")
#
# set.seed(2016)
# ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ] # small data for fast test
# set.seed(2016)
# # # TODO bug in gbm: following causes error: survived is not of type numeric, ordered, or factor
# # ptit$survived <- ptit$survived == "survived"
# ptit <- ptit[!is.na(ptit$age), ]
# # TODO change this to build same model as gbm.gaussian
# train_params <-
#      training_params(num_trees = 50,
#                      shrinkage = 0.1,
#                      bag_fraction = 0.5,
#                      num_train = round(.8 * nrow(ptit)))
# par(mfrow=c(2,2), mar=c(3,3,4,1))
# set.seed(2016)
# gbmt.gaussian <- gbmt(age~., data=ptit,
#             distribution=gbm_dist("Gaussian"),
#             train_params = train_params,
#             is_verbose = FALSE)
# expect.err(try(plotres(gbmt.gaussian)),
#            "use keep.data=TRUE in the call to gbm")
# set.seed(2016)
# gbmt.gaussian <- gbmt(age~., data=ptit,
#             distribution=gbm_dist("Gaussian"),
#             train_params = train_params,
#             is_verbose = FALSE,  keep_gbm_data=TRUE)
# w1 <- plotres(gbmt.gaussian, which=1, do.par=FALSE, w1.smooth=TRUE,
#               w1.main="gbmt.gaussian")
# cat("w1 plot for gbmt.gaussian returned (w1.smooth=TRUE):\n")
# print(w1)
# plot(0, 0) # dummy plot
# set.seed(2016)
# w3 <- plotres(gbmt.gaussian, which=3, do.par=FALSE, info=TRUE,
#         smooth.col=0, col=ptit$sex, # ylim=c(-40,40),
#         wmain="nresponse=1")
#
# # compare to manual residuals
# iused <- 1:(train.frac * nrow(ptit))
# y <- ptit$age[iused]
# n.trees <- plotmo:::gbm.n.trees(gbmt.gaussian)
# # TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# # yhat <- predict(gbmt.gaussian, type="response", n.trees=n.trees)
# yhat <- predict(gbmt.gaussian, newdata=ptit, type="response", n.trees=n.trees)
# yhat <- yhat[iused]
# plot(yhat, y - yhat,
#      col=ptit$sex[iused], main="manual gaussian residuals",
#      pch=20, ylim=c(-40,40))
# abline(h=0, col="gray")
# stopifnot(all(yhat == w3$x))
# stopifnot(all(y - yhat == w3$y))
# par(org.par)
#
# w1 <- plotres(gbmt.gaussian, predict.n.trees=13, w1.grid.col=1, trace=1, SHOWCALL=TRUE,
#               w1.smooth=TRUE,
#               w1.main="predict.n.trees=13 w1.grid.col=1")
# cat("second w1 plot for gbmt.gaussian returned (w1.smooth=TRUE):\n")
# print(w1)
# plotmo(gbmt.gaussian, trace=-1, SHOWCALL=TRUE)
#
# par(org.par)
#
# cat("--- distribution=\"bernoulli\" ----------------------------------\n")
#
# set.seed(2016)
# ptit <- ptitanic[sample(1:nrow(ptitanic), size=80), ]
# ptit$survived <- ptit$survived == "survived"
# temp <- ptit$pclass # put pclass at the end so can check ordering of importances
# ptit$pclass <- NULL
# ptit$pclass <- factor(as.numeric(temp), labels=c("first", "second", "third"))
# # TODO change this to build same model as gbm.bernoulli
# train_params <-
#      training_params(num_trees = 100,
#                      shrinkage = 0.1,
#                      bag_fraction = 0.5,
#                      num_train = round(.8 * nrow(ptit)))
# set.seed(2016)
# gbmt.bernoulli <- gbmt(survived~., data=ptit,
#             distribution=gbm_dist("Bernoulli"),
#             train_params = train_params,
#             cv_folds = 3,
#             is_verbose = FALSE,  keep_gbm_data=TRUE)
# par(mfrow=c(2,2))
# par(mar=c(3.5, 3, 2, 0.5))  # small margins and text to pack figs in
# par(mgp=c(1.5, .4, 0))      # squash axis annotations
# w1 <- plotres(gbmt.bernoulli, which=c(1,4),
#         col=ptit$survived+2, trace=0, do.par=FALSE,
#         w1.main="gbmt.bernoulli")
# cat("w1 plot for gbmt.bernoulli with cv.folds=3 returned:\n")
# print(w1)
#
# w3 <- plotres(gbmt.bernoulli, which=3, predict.n.trees=40,
#         ylim=c(-.6, 1), xlim=c(.1, .6),
#         col=ptit$sex, trace=0, do.par=FALSE, smooth.col=0)
#
# # compare to manual residuals
# iused <- 1:(train.frac * nrow(ptit))
# y <- ptit$survived[iused]
# # TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# # yhat <- predict(gbmt.bernoulli, type="response", n.trees=40)
# yhat <- predict(gbmt.bernoulli, newdata=ptit, type="response", n.trees=40)
# yhat <- yhat[iused]
# plot(yhat, y - yhat, col=ptit$sex,
#      main="manual bernoulli residuals", pch=20, cex=1,
#      ylim=c(-.6, 1), xlim=c(.1, .6))
# abline(h=0, col="gray")
# stopifnot(all(yhat == w3$x))
# stopifnot(all(y - yhat == w3$y))
# par(org.par)
#
# plotmo(gbmt.bernoulli, do.par=2)
# print(summary(gbmt.bernoulli)) # will also plot
# par(org.par)

cat("--- gbm3: distribution=\"gaussian\", formula interface ----------------------------------\n")

library(gbm3)

set.seed(2016)
ptit <- ptitanic[sample(1:nrow(ptitanic), size=70), ] # small data for fast test
set.seed(2016)
# # TODO bug in gbm3: following causes error: survived is not of type numeric, ordered, or factor
# ptit$survived <- ptit$survived == "survived"
ptit <- ptit[!is.na(ptit$age), ]
train.frac <- .8
set.seed(2016)
gbm3.gaussian <- gbm3::gbm(age~., data=ptit, train.frac=train.frac,
                   distribution="gaussian",
                   n.trees=50, shrinkage=.1, keep.data=FALSE)
expect.err(try(plotres(gbm3.gaussian)), "use keep_gbm_data=TRUE in the call to gbm")
set.seed(2016)
gbm3.gaussian <- gbm3::gbm(age~., data=ptit, train.frac=train.frac,
                   distribution="gaussian",
                   n.trees=50, shrinkage=.1)
par(mfrow=c(2,2), mar=c(3,3,4,1))
w1 <- plotres(gbm3.gaussian, which=1, do.par=FALSE, w1.smooth=TRUE,
        w1.main="gbm3.gaussian")
cat("w1 plot for gbm3.gaussian returned (w1.smooth=TRUE):\n")
print(w1)
plot(0, 0) # dummy plot
w3 <- plotres(gbm3.gaussian, which=3, do.par=FALSE, info=TRUE,
        smooth.col=0, col=ptit$sex, # ylim=c(-40,40),
        wmain="nresponse=1")

# compare to manual residuals
iused <- 1:(train.frac * nrow(ptit))
y <- ptit$age[iused]
n.trees <- plotmo:::gbm.n.trees(gbm3.gaussian)
# TODO following fails in the new version of gbm (version 2.2) (you have to provide newdata)
# yhat <- predict(gbm3.gaussian, type="response", n.trees=n.trees)
yhat <- predict(gbm3.gaussian, newdata=ptit, type="response", n.trees=n.trees)
yhat <- yhat[iused]
plot(yhat, y - yhat,
     col=ptit$sex[iused], main="manual gaussian residuals",
     pch=20, ylim=c(-40,40))
abline(h=0, col="gray")
stopifnot(all(yhat == w3$x))
stopifnot(all(y - yhat == w3$y))
par(org.par)

w1 <- plotres(gbm3.gaussian, predict.n.trees=13, w1.grid.col=1, trace=1, SHOWCALL=TRUE,
              w1.smooth=TRUE,
              w1.main="predict.n.trees=13 w1.grid.col=1")
cat("second w1 plot for gbm3.gaussian returned (w1.smooth=TRUE):\n")
print(w1)
plotmo(gbm3.gaussian, trace=-1, SHOWCALL=TRUE)
# plotmo(gbm3.gaussian, trace=-1, all1=TRUE, SHOWCALL=TRUE)
# plotmo(gbm3.gaussian, trace=-1, all2=TRUE, SHOWCALL=TRUE)

cat("--- gbm3: distribution=\"gaussian\", xy interface ----------------------------------\n")

y = ptit$age
x = ptit[,c(1,2,3,5,6)]
train_params=gbm3::training_params(num_trees=100,
 interaction_depth=2,
 min_num_obs_in_node=3,
 shrinkage=0.1, bag_fraction=0.5,
 id=seq_len(nrow(x)), num_train=round(0.5 * nrow(x)),
 num_features=ncol(x))
gbm3fit <- gbm3::gbmt_fit(x, y, train_params=train_params,
                          keep_gbm_data=TRUE, dist=gbm_dist("Gaussian"))
plotmo(gbm3fit, trace=-1, SHOWCALL=TRUE)
plotres(gbm3fit, trace=-1, SHOWCALL=TRUE)

cat("--- gbm3: large number of variables ----------------------------------\n")

set.seed(2024)
N <- 1000

X <- data.frame(X1=runif(N), X2=2*runif(N), X3=3*runif(N),
                X4=runif(N), X5=2*runif(N), X6=3*runif(N),
                X7=runif(N), X8=2*runif(N), X9=3*runif(N),
                X10=runif(N), X11=2*runif(N), X12=3*runif(N),
                X13=runif(N), X14=2*runif(N), X15=3*runif(N))

# Y <- sample(c(0, 1), N, replace = TRUE)
set.seed(2024)
Y <- sqrt(X[,1])  +
     sqrt(X[,2])  +
     sqrt(X[,3])  +
     sqrt(X[,4])  +
     sqrt(X[,5])  +
     sqrt(X[,6])  +
     .5 * sqrt(X[,8]) +
     sqrt(X[,9])  +
     sqrt(X[,10]) +
     sqrt(X[,11]) +
     sqrt(X[,12])

data <- data.frame(Y, X)
set.seed(2024)
gbm3.big <- gbm3::gbm(Y~., data=data, shrinkage=0.1, dist="gaussian")
y = data[,1]
x = data[,2:ncol(data)]
train_params=gbm3::training_params(num_trees=100,
 interaction_depth=3,
 min_num_obs_in_node=10,
 shrinkage=0.1, bag_fraction=0.5,
 id=seq_len(nrow(x)), num_train=round(0.5 * nrow(x)),
 num_features=ncol(x))
gbm3fit.big <- gbm3::gbmt_fit(x, y, train_params=train_params, keep_gbm_data=TRUE, dist=gbm_dist("Gaussian"))

set.seed(2024)
plotmo(gbm3.big, SHOWCALL=TRUE)
plotmo(gbm3.big, all1=TRUE, all2=TRUE, caption="all1=TRUE, all2=TRUE")
plotmo(gbm3.big, all1=TRUE, all2=2, caption="all1=TRUE, all2=2")
plotres(gbm3.big, trace=-1, SHOWCALL=TRUE)

set.seed(2024)
plotmo(gbm3fit.big, SHOWCALL=TRUE)
plotmo(gbm3fit.big, all1=TRUE, caption="all1=TRUE")
plotmo(gbm3fit.big, all2=TRUE, caption="all2=TRUE")
plotmo(gbm3fit.big, all2=2, caption="all2=2")
plotres(gbm3.big, trace=-1, SHOWCALL=TRUE)

source("test.epilog.R")
