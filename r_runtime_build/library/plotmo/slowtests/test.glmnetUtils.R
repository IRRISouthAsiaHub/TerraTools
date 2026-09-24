# test.glmnet.R: glmnetUtils tests for plotmo and plotres

source("test.prolog.R")
library(earth)
library(glmnetUtils)
data(ozone1)
data(etitanic)
get.tit <- function() # abbreviated titanic data
{
    tit <- etitanic
    pclass <- as.character(tit$pclass)
    # change the order of the factors so not alphabetical
    pclass[pclass == "1st"] <- "first"
    pclass[pclass == "2nd"] <- "class2"
    pclass[pclass == "3rd"] <- "classthird"
    tit$pclass <- factor(pclass, levels=c("class2", "classthird", "first"))
    # log age is so we have a continuous predictor even when model is age~.
    set.seed(2015)
    tit$logage <- log(tit$age) + rnorm(nrow(tit))
    tit$parch <- NULL
    # by=12 gives us a small fast model with an additive and a interaction term
    tit <- tit[seq(1, nrow(etitanic), by=12), ]
}
plotmores <- function(object, ..., trace=0, SHOWCALL=TRUE, title.extra="", ncol=2) {
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(2,ncol))
    caption <- paste(deparse(substitute(object)), collapse=" ")
    call <- match.call(expand.dots=TRUE)
    call <- strip.space(paste(deparse(substitute(call)), collapse=" "))
    call <- gsub(",", ", ", call)
    call <- paste(title.extra, call, sep="")
    printf("%s\n", call)
    # plotmo on glmnet mods is boring but we test it anyway
    plotres(object, trace=trace, SHOWCALL=SHOWCALL, do.par=FALSE, which=c(1,3), ...)
    title(paste("\n", call), outer=TRUE)
    plotmo(object, trace=trace, SHOWCALL=SHOWCALL, do.par=FALSE, ...)
}
tit <- get.tit()
set.seed(2015)
xmat <- as.matrix(tit[,c(2,5,6)])
agedata <- data.frame(tit[,4], xmat)
colnames(agedata) <- c("age", "survived", "sibsp", "logage")
set.seed(2015)
mod.glmnet.xmat <- glmnet(xmat, tit[,4]) # tit[,4] is age
plotres(mod.glmnet.xmat)
plotmo(mod.glmnet.xmat)
plotmores(mod.glmnet.xmat, predict.s=2.5)

mod.glmnet.agedata <- glmnet(age~., data=agedata)
expect.err(try(plotres(mod.glmnet.agedata)), "for this plot, glmnet.formula must be called with use.model.frame=TRUE")
mod.glmnet.agedata <- glmnet(age~., data=agedata, use.model.frame=TRUE)
plotmores(mod.glmnet.agedata, predict.s=2.5)

set.seed(2015)
mod.cv.glmnet.xmat <- cv.glmnet(xmat, tit[,4], nfolds=3)

cat("==Test plotmo trace=1 and lambda.min\n")
plotmores(mod.cv.glmnet.xmat, predict.s="lambda.min", trace=1, ncol=3)

set.seed(2015)
mod.cv.glmnet.agedata <- cv.glmnet(age~., data=agedata)
expect.err(try(plotres(mod.cv.glmnet.agedata)), "for this plot, cv.glmnet.formula must be called with use.model.frame=TRUE")
set.seed(2015)
mod.cv.glmnet.agedata <- cv.glmnet(age~., data=agedata, use.model.frame=TRUE)
cat("==Test lambda.min\n")
plotmores(mod.cv.glmnet.agedata, predict.s="lambda.min", trace=1, ncol=3)

printf("======== binomial model\n")

set.seed(2016)
n <- 50
p <- 4
xx <- matrix(rnorm(n*p), n, p)
colnames(xx) <- paste("x", 1:ncol(xx), sep="")
yy <- ifelse(xx[,1] + xx[,2] + rnorm(n) > .5, TRUE, FALSE)
print(cov(xx, yy))
yy <- factor(yy)
dataxy <- data.frame(yy, xx)
binomial.mod <- glmnet(xx, yy, family="binomial")
plotmores(binomial.mod, ncol=3)
binomial.mod.form <- glmnet(yy~., data=dataxy, family="binomial", use.model.frame=TRUE)
plotmores(binomial.mod.form, ncol=3)
par(org.par)

printf("======== glmnet family=\"mgaussian\"\n")
set.seed(2015)
p <- 10
n <- 30
xx <- cbind((1:n)/n, matrix(rnorm(n*(p-1)),n,p-1))
colnames(xx) <- paste0("x", 1:p)
# ymultresp <- cbind(rowSums(xx[,1:5]^3), rowSums(xx[,5:p]^3), 1:n)
set.seed(1)
ymultresp <- cbind(xx[,1]+.001*rnorm(n), rowSums(xx[,2:5]^3), rnorm(n))
glmnet.mgaussian <- glmnet(xx, ymultresp, family="mgaussian")
plotres(glmnet.mgaussian, nresponse=1, SHOWCALL=TRUE, which=c(1:3), do.par=2, info=1)
# manually calculate the residuals
plot(x=predict(glmnet.mgaussian, newx=xx, s=0)[,1,1],
     y=ymultresp[,1] - predict(glmnet.mgaussian, newx=xx, s=0)[,1,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manually calculated residuals, nresponse=1, s=0")
abline(h=0, col="gray")
par(org.par)

# # TODO is glmnet mgaussian supported with a formula interface?
# dataxy <- data.frame(ymultresp, xx)
# colnames(dataxy) <- c("y1", "y2", "y3", "x1", "x2", "x3", "x4", "x5", "x5", "x6", "x7", "x8", "x9", "x10")
# glmnet.mgaussian.form <- glmnet(xx, ymultresp, family="mgaussian")
# plotres(glmnet.mgaussian.form, nresponse=1, SHOWCALL=TRUE, which=c(1:3), do.par=2, info=1)

par(mfrow=c(2,3), mar=c(3,3,3,.5), oma=c(0,0,3,0), mgp=c(1.5,0.4,0), tcl=-0.3)

data(trees)
set.seed(2015)
# variable with a long name
x50 <- cbind(trees[,1:2], Girth12345678901234567890=rnorm(nrow(trees)))
mod.with.long.name <- glmnet(data.matrix(x50),data.matrix(trees$Volume))
plotmores(mod.with.long.name, ncol=3)
data.x50 <- data.frame(trees$Volume, x50)
colnames(data.x50) <- c("Volume",  "Girth", "Height", "Girth12345678901234567890")
mod.with.long.name.form <- glmnet(Volume~., data=data.x50, use.model.frame=TRUE)
plotmores(mod.with.long.name.form, ncol=3)
par(org.par)

#-- make sure that we can work with all families

set.seed(2016)
par(mfrow=c(3,3), mar=c(3,3,3,1))
n <- 100
p <- 4
xx <- matrix(rnorm(n*p), n, p)
g2 <- sample(1:2, n, replace=TRUE)
data.xg2 <- data.frame(g2, xx)
for(family in c("gaussian","binomial","poisson")) {
    title.extra <- paste(family, ": ")
    mod <- glmnet(xx,g2,family=family)
    plotmores(mod, xvar="lambda", ncol=3, title.extra=title.extra)
    title.extra <- paste("formula", family, ": ")
    mod.form <- glmnet(g2~., data.xg2, family=family, use.model.frame=TRUE)
    plotmores(mod.form, xvar="lambda", ncol=3, title.extra=title.extra)
}
par(org.par)
# cox
library(plotmo)
n <- 100
p <- 20
nzc <- trunc(p/10)
set.seed(2016)
beta <- rnorm(nzc)
x7 <- matrix(rnorm(n*p), n, p)
beta <- rnorm(nzc)
fx <- x7[,seq(nzc)] %*% beta/3
hx <- exp(fx)
ty <- rexp(n, hx)
tcens <- rbinom(n=n, prob=.3, size=1)# censoring indicator
yy <- cbind(time=ty, status=1-tcens) # yy=Surv(ty,1-tcens) with library(survival)
glmnet.cox <- glmnet(x=x7, y=yy, family="cox")
plotmores(glmnet.cox, ncol=3, degree1=1:4)
par(org.par)
# TODO formula interface not tested for cox models

source("test.epilog.R")
