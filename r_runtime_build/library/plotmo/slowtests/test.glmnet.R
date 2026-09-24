# test.glmnet.R: glmnet tests for plotmo and plotres

source("test.prolog.R")
library(earth)
library(glmnet)
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
tit <- get.tit()
set.seed(2015)
xmat <- as.matrix(tit[,c(2,5,6)])
set.seed(2015)
mod.glmnet.xmat <- glmnet(xmat, tit[,4])
# plotmo on glmnet mods is boring but we test it anyway
plotmo1(mod.glmnet.xmat)
plotres1(mod.glmnet.xmat)

# compare to plot.glmnet
par(mfrow=c(4,2), mar=c(3,6,3.5,6)) # extra side margins for more square plots
plot_glmnet(mod.glmnet.xmat, main="mod.glmnet.xmat\ncompare to plot.glmnet")
plot(0,0)
plot_glmnet(mod.glmnet.xmat, xvar="norm", col=c(3,2,1))
plot(mod.glmnet.xmat, xvar="norm")
plot_glmnet(mod.glmnet.xmat, xvar="lambda")
plot(mod.glmnet.xmat, xvar="lambda")
plot_glmnet(mod.glmnet.xmat, xvar="dev")
plot(mod.glmnet.xmat, xvar="dev")
par(org.par)

set.seed(2015)
mod.cv.glmnet.xmat <- cv.glmnet(xmat, tit[,4], nfolds=3)

# following was needed before plotmo 3.1.3 (before adding plotmo.prolog.cv.glmnet)
# mod.cv.glmnet.xmat$x <- as.data.frame(xmat)
# mod.cv.glmnet.xmat$y <- tit[,4]

cat("==Test plotmo trace=1 and lambda.min\n")
plotmo1(mod.cv.glmnet.xmat, predict.s="lambda.min", trace=1)
cat("==Test plotmo trace=2 and lambda.min\n")
plotmo1(mod.cv.glmnet.xmat, predict.s="lambda.min", trace=2)
cat("==Test plotres trace=1 and lambda.1se\n")
plotres1(mod.cv.glmnet.xmat, predict.s="lambda.1se", trace=1)
cat("==Test plotres trace=2 and lambda.1se\n")
plotres1(mod.cv.glmnet.xmat, predict.s="lambda.1se", trace=2)

set.seed(2015)
x <- matrix(rnorm(100*20),100,20) # 20 variables
y <- rnorm(100)
mod <- glmnet(x,y)
plotmo1(mod)

# test w1.label
par(mfrow=c(2,3))
par(cex=1)
par(mar=c(3,3,3,1))
plotres(mod, which=1,                w1.main="default w1.label")
plotres(mod, which=1, w1.label=5,    w1.main="w1.label=5")
plotres(mod, which=1, w1.label=0,    w1.main="w1.label=0")
plotres(mod, which=1, w1.label=TRUE, w1.main="w1.label=TRUE")
plotres(mod, which=1, w1.label=100,  w1.main="w1.label=100")
par(org.par)

# test w1 and non w1 args passed
par(mfrow=c(2,2), mar=c(4,4,4,4), cex=1)

plot_glmnet(mod, w1.col=3:4, w1.xvar="norm",
            main="plot_glmnet\nw1.col=3:4 w1.xvar=\"norm\"")

plot_glmnet(mod, col=3:4, xvar="norm",
            main="plot_glmnet\ncol=3:4 xvar=\"norm\"")

plot_glmnet(mod, col=3:4, w1.col=1:2,
            w1.xvar="norm", xvar="lambda",
            main="plot_glmnet\ncol=3:4 w1.col=1:2\nw1.xvar=\"norm\", xvar=\"lambda\"")

par(org.par)
par(mfrow=c(3,2), mar=c(3,4,4,4), cex=1)

plotres(mod, which=c(1,3), do.par=FALSE, w1.col=3:4, w1.xvar="norm",
        w1.main="plotres\nw1.col=3:4 w1.xvar=\"norm\"")

plotres(mod, which=c(1,3), do.par=FALSE, col=3:4, xvar="norm",
        w1.main="plotres\nplotres\ncol=3:4 xvar=\"norm\"")

plotres(mod, which=c(1,3), do.par=FALSE, col=3:4, w1.col=1:2,
        w1.main="plotres\ncol=3:4 w1.col=1:2")

par(org.par)

# glmnet with sparse matrices
set.seed(2015)
n <- 100
p <- 20
nzc <- trunc(p/10)
x <- matrix(rnorm(n*p),n,p)
iz <- sample(1:(n*p),size=n*p*.85,replace=FALSE)
x[iz] <- 0
sx <- Matrix(x,sparse=TRUE)
# colnames(sx) <- paste("x", 1:ncol(sx), sep="") # need column names for plotmo
inherits(sx,"sparseMatrix") # confirm that it is sparse
beta <- rnorm(nzc)
fx <- x[,seq(nzc)]%*%beta
eps <- rnorm(n)
y <- fx+eps
px <- exp(fx)
px <- px/(1+px)
ly <- rbinom(n=length(px),prob=px,size=1)
mod.glmnet.sx <- glmnet(sx,y)
plotmo1(mod.glmnet.sx, all2=TRUE) # will give warning: too many predictors to plot all pairs
plotmo1(mod.glmnet.sx, all2=2, caption="all2=2") # test all2=2
plotmo1(mod.glmnet.sx, all2=2, degree2=1:3, caption="all2=2 degree2=1:3")
plotres(mod.glmnet.sx)
par(org.par)

par(mfrow=c(2,4), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
y <- trees$Volume
x <- as.matrix(data.frame(Girth=trees$Girth, Height=trees$Height))
glmnet <- glmnet(x, y)
plotres(glmnet, do.par=FALSE, caption="glmnet and lm: top and bottom should be the same")
lm <- lm(Volume~., data=trees)
plotres(lm, do.par=FALSE, SHOWCALL=TRUE)

par(mfrow=c(3,2), mar=c(3,3,3,1), mgp=c(1.5,0.5,0), oma=c(0,0,2.5,0))
plotres(glmnet, do.par=FALSE, which=c(1,3), w1.xvar="norm",
        caption="glmnet with various options", SHOWCALL=TRUE)
plotres(glmnet, trace=1, do.par=FALSE, which=c(1,3), SHOWCALL=TRUE)
plotres(glmnet, trace=1, do.par=FALSE, which=c(1,3), predict.s=5, SHOWCALL=TRUE)
par(org.par)

printf("======== glmnet additional tests\n")
set.seed(2015)
p <- 10
n <- 30
x <- cbind(matrix(rnorm(n*p),n,p))
y <- rowSums(x[,1:3]^3)
glmnet <- glmnet(x,y)
plotres(glmnet, SHOWCALL=TRUE, caption="glmnet: y <- rowSums(x[,1:3]^3)")
plotres(glmnet, SHOWCALL=TRUE, w1.xvar="norm")
par(mfrow=c(1,1))
omar <- par("mar")
ocex.axis <- par("cex.axis")
ocex.lab <- par("cex.lab")
plotres(glmnet, SHOWCALL=TRUE, which=1)
stopifnot(par("mar") == omar)
stopifnot(par("cex.axis") == ocex.axis)
stopifnot(par("cex.lab") == ocex.lab)
par(org.par)

# test some args for plot_glmnet
plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, trace=0, col.main=2,
        w1.xlab="my xlab", w1.ylab="my ylab",
        w1.main="test some args for plot_glmnet1",
        w1.col=4:1)

plot_glmnet(glmnet, trace=0, col.main=2, main="test some args for plot_glmnet2",
        xlab="my xlab", ylab="my ylab",
        col=4:1, ylim=c(-2,4)) # TODO xlim=c(-5,3))

plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, which=c(1,3), grid.col="gray", do.par=2)
plotres(glmnet, predict.s=.05, SHOWCALL=TRUE, which=c(1,3), w1.s.col=0, do.par=0)
par(org.par)

# TODO the following issues a stream of warnings: restarting interrupted promise evaluation
expect.err(try(plotres(glmnet, w1.col=nonesuch)), "cannot evaluate 'col'")

printf("======== glmnet multinomial (multnet)\n")
par(mfrow=c(4,4), mar=c(3,3,3,1))
set.seed(2016)
n <- 200
p <- 4
x <- matrix(rnorm(n*p), n, p)
colnames(x) <- paste("x", 1:ncol(x), sep="")

# "1" is correlated with x[,1], "4" is correlated with x[,2], "2" and "3" not correlated
y <- ifelse(x[,1] > 0.5, 1,
     ifelse(x[,2] > 0.0, 4,
     sample(c(2,3), size=nrow(x), replace=TRUE)))
print(cov(x, y))
y <- factor(y)

# TODO Following causes the following warning:
#      Warning: from glmnet Fortran code (error code -90); Convergence for 90th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
multinomial.mod <- glmnet(x, y, family="multinomial")

plotres(multinomial.mod, nresponse=1, w1.main="nresponse=1",
        main="family=\"multinomial\"",
        smooth.col=0, info=TRUE,
        trace=0, which=c(1,3), do.par=FALSE, xlim=c(-.2,1.2), ylim=c(-1.2, 1.2))

plotres(multinomial.mod, nresponse=2, w1.main="nresponse=2",
        smooth.col=0, info=TRUE,
        trace=0, which=c(1,3), do.par=FALSE, xlim=c(-.2,1.2), ylim=c(-1.2, 1.2))

plotres(multinomial.mod, nresponse=3, w1.main="nresponse=3",
        smooth.col=0, info=TRUE,
        trace=0, which=c(1,3), do.par=FALSE, xlim=c(-.2,1.2), ylim=c(-1.2, 1.2))

plotres(multinomial.mod, nresponse=4, w1.main="nresponse=4",
        smooth.col=0, info=TRUE,
        trace=0, which=c(1,3), do.par=FALSE, xlim=c(-.2,1.2), ylim=c(-1.2, 1.2))

plotmo(multinomial.mod, nresponse=1, trace=0, do.par=FALSE, degree1=1:2)
plotmo(multinomial.mod, nresponse=2, trace=0, do.par=FALSE, degree1=1:2)

par(mgp=c(1.5, .4, 0))
plot(multinomial.mod, xvar="norm") # compare to plot.glmnet
par(org.par)

# compare to earth
par(mfrow=c(4,3), mar=c(3,3,1,1))
yfac <- factor(c("a","b","c","d")[y])
earth.mod <- earth(x, yfac, trace=0)

plotres(earth.mod, nresponse=1,
        main=sprint("multiresponse\nnresponse=1   rsq %.2g", earth.mod$rsq.per.response[1]),
        which=3, xlim=c(-.2, 1.2), ylim=c(-1.2, 1.2),
        smooth.col=0, info=TRUE,
        do.par=FALSE, trace=0, jitter=7, cex.response=.7)
plotmo(earth.mod, nresponse=1, do.par=FALSE)

plotres(earth.mod, nresponse=2,
        main=sprint("nresponse=2   rsq %.2g", earth.mod$rsq.per.response[2]),
        which=3, xlim=c(-.2, 1.2), ylim=c(-1.2, 1.2),
        smooth.col=0, info=TRUE,
        do.par=FALSE, trace=0, jitter=7, cex.response=.7)
plotmo(earth.mod, nresponse=2, do.par=FALSE)

plotres(earth.mod, nresponse=3,
        main=sprint("nresponse=3   rsq %.2g", earth.mod$rsq.per.response[3]),
        which=3, xlim=c(-.2, 1.2), ylim=c(-1.2, 1.2),
        smooth.col=0, info=TRUE,
        do.par=FALSE, trace=0, jitter=7, cex.response=.7)
plotmo(earth.mod, nresponse=3, do.par=FALSE)

plotres(earth.mod, nresponse=4,
        main=sprint("nresponse=4   rsq %.2g", earth.mod$rsq.per.response[4]),
        which=3, xlim=c(-.2, 1.2), ylim=c(-1.2, 1.2),
        smooth.col=0, info=TRUE,
        do.par=FALSE, trace=0, jitter=7, cex.response=.7)
plotmo(earth.mod, nresponse=4, do.par=FALSE)

print(summary(earth.mod))

par(org.par)

printf("======== binomial model\n")

set.seed(2019)
n <- 50
p <- 4
x <- matrix(rnorm(n*p), n, p)
colnames(x) <- paste("x", 1:ncol(x), sep="")
y <- ifelse(x[,1] + x[,2] + .1 * rnorm(n) > .5, TRUE, FALSE)
print(cov(x, y))
y <- factor(y)
glmnet.binomial <- glmnet(x, y, family="binomial")
par(mfrow=c(2,3), mar=c(3,3,1,1))
plotres(glmnet.binomial, info=T, predict.s=.02, which=c(1,3), do.par=FALSE, w1.main="glmnet.binomial")
plot(glmnet.binomial)
earth.mod <- earth(x, y)
set.seed(2019)
plotres(earth.mod, info=T, which=c(1,3), do.par=FALSE)
par(org.par)
par(mfrow=c(2,4), mar=c(3,3,1,1))
set.seed(2019)
plotmo(glmnet.binomial, do.par=FALSE)
plotmo(earth.mod, do.par=FALSE, main="binomial earth.mod")
par(org.par)

printf("======== glmnet family=\"mgaussian\"\n")
set.seed(2015)
p <- 10
n <- 30
x <- cbind((1:n)/n, matrix(rnorm(n*(p-1)),n,p-1))
colnames(x) <- paste0("x", 1:p)
# ymultresp <- cbind(rowSums(x[,1:5]^3), rowSums(x[,5:p]^3), 1:n)
set.seed(1)
ymultresp <- cbind(x[,1]+.001*rnorm(n), rowSums(x[,2:5]^3), rnorm(n))
glmnet.mgaussian <- glmnet(x, ymultresp, family="mgaussian")
plotres(glmnet.mgaussian, nresponse=1, SHOWCALL=TRUE, which=c(1:3), do.par=2, info=1)
# manually calculate the residuals
plot(x=predict(glmnet.mgaussian, newx=x, s=0)[,1,1],
     y=ymultresp[,1] - predict(glmnet.mgaussian, newx=x, s=0)[,1,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manually calculated residuals, nresponse=1, s=0")
abline(h=0, col="gray")
par(org.par)
plotres(glmnet.mgaussian, nresponse=2, SHOWCALL=TRUE, which=c(1:3), do.par=2, info=1)
# manually calculate the residuals
plot(x=predict(glmnet.mgaussian, newx=x, s=0)[,2,1],
     y=ymultresp[,2] - predict(glmnet.mgaussian, newx=x, s=0)[,2,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manually calculated residuals, nresponse=2, s=0")
abline(h=0, col="gray")
par(org.par)
plotmo(glmnet.mgaussian, nresponse=1, SHOWCALL=TRUE)
plotmo(glmnet.mgaussian, nresponse=2, SHOWCALL=TRUE)

graphics::par(mfrow=c(2,2), mgp=c(1.5,0.4,0), tcl=-0.3, cex.main=1,
              font.main=1, mar=c(4,3,1.2,0.8), oma=c(0,0,4,0), cex=0.83)

plotres(glmnet.mgaussian, nresponse=2, SHOWCALL=TRUE, which=3, do.par=FALSE,
        caption="glmnet.mgaussian compare to manually calculated residuals")
plot(x=predict(glmnet.mgaussian, newx=x, s=0)[,2,1],
     y=ymultresp[,2] - predict(glmnet.mgaussian, newx=x, s=0)[,2,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manual residuals, nresponse=2, s=0")
abline(h=0, col="gray")

plotres(glmnet.mgaussian, nresponse=2, predict.s=.5, SHOWCALL=TRUE, which=3, do.par=FALSE)
plot(x=predict(glmnet.mgaussian, newx=x, s=.5)[,2,1],
     y=ymultresp[,2] - predict(glmnet.mgaussian, newx=x, s=.5)[,2,1],
     pch=20, xlab="Fitted", ylab="Residuals",
     main="Manual residuals, nresponse=2, s=.5")
abline(h=0, col="gray")

plotres(glmnet.mgaussian, predict.s=.05, nresponse=3, info=TRUE, SHOWCALL=TRUE) # essentially random

par(org.par)
par(mfrow=c(2,3), mar=c(3,3,3,.5), oma=c(0,0,3,0), mgp=c(1.5,0.4,0), tcl=-0.3)

data(trees)
set.seed(2015)
# variable with a long name
x50 <- cbind(trees[,1:2], Girth12345678901234567890=rnorm(nrow(trees)))
mod.with.long.name <- glmnet(data.matrix(x50),data.matrix(trees$Volume))
plotres(mod.with.long.name, which=1, caption="test plot_glmnet with x50 and x60")

# one inactive variable (all coefs are zero for variable "rand")
set.seed(2015)
x60 <- cbind(trees[,1], rand=rnorm(nrow(trees)), trees[,2])
# complicate the issue: use an unnamed column (column 3)
colnames(x60) <- c("Girth", "rand", "")
mod.with.inactive.var <- glmnet(data.matrix(x60),data.matrix(trees$Volume))
mod.with.inactive.var$beta["rand",] = 0 # TODO hack force inactive variable
plotres(mod.with.inactive.var, which=1)
plotres(mod.with.inactive.var, which=1, w1.xvar="norm")
# compare to plot.glmnet (but note that labels aren't always plotted unless par=c(1,1)?)
plot(mod.with.inactive.var, xvar="norm", label=TRUE)
# plotmo calls the unnamed column "x3", fair enough
plotmo(mod.with.inactive.var, do.par=FALSE, pt.col=2)

# single active variable
x70 <- cbind(trees[,1,drop=F], 0)
a <- glmnet(data.matrix(x70), data.matrix(trees$Volume))
par(org.par)
par(mfrow=c(2,2), mar=c(3,3,2,4))
plotres(a, which=1, predict.s=1, caption="single active variable")
plotres(a, which=1, w1.xvar="norm")
plotres(a, which=1, w1.xvar="lambda")
plotres(a, which=1, w1.xvar="dev")

#--- test interaction of w1. and non w1 args -------------------------------------

#--- glmnet model, which=1 ---

par(org.par)
par(mfrow=c(4,3), mar=c(3, 3, 4, 1), mgp=c(2, 0.6, 0))

plotres(mod.glmnet.xmat, which=1,
        w1.xlim=c(6,-6),
        w1.ylim=c(-5,5),
        w1.col=1:2,
        w1.main="TEST INTERACTION OF W1 ARGS PAGE 1 (which=1)\n\nwhich=1 w1.xlim=c(6,-6)\nw1.ylim=c(-5,5)) w1.col=1:2,")

plotres(mod.glmnet.xmat, which=1, cex.main=1.2,
        xlim=c(9,-9),
        ylim=c(-60,60),
        col=3:4,
        w1.main="which=1 xlim=c(9,-9)\nylim=c(-60,60)) col=3:4,")

plotres(mod.glmnet.xmat, which=1, cex.main=1,
        xlim=c(9,-9), w1.xlim=c(6,-6),
        ylim=c(-60,60), w1.ylim=c(-5,5),
        w1.col=1:2, col=3:4,
        w1.main="which=1 xlim=c(9,-9), w1.xlim=c(6,-6)\nylim=c(-60,60), w1.ylim=c(-5,5)) w1.col=1:2, col=3:4")

#--- glmnet model, which=c(1,3,4) ---

plotres(mod.glmnet.xmat, which=c(1,3,4), cex.main=1,
        ylim=c(-70,70),    xlim=c(-20, 60),
        col=2:3, do.par=FALSE,
        w1.main="TEST INTERACTION OF W1 ARGS PAGE 1 (which=c(1,3,4))\nlim=c(-70,70), xlim=c(-20, 60)")

plotres(mod.glmnet.xmat, which=c(1,3,4), cex.main=1.2,
        ylim=c(-70,70),    xlim=c(-20, 60), qq.xlim=c(-7,5),
        col=2:3, do.par=FALSE,
        w1.main="ylim=c(-70,70), xlim=c(-20, 60)\nqq.xlim=c(-7,5)")

plotres(mod.glmnet.xmat, which=c(1,3,4), cex.main=1.2,
        w1.ylim=c(-7,7), w1.xlim=c(4,-4),     col=2:3, do.par=FALSE,
        w1.main="w1.ylim=c(-7,7), w1.xlim=c(4,-4)")

# plotres(mod.glmnet.xmat, which=c(1,3,4), cex.main=.9,
#         w1.ylim=c(-7,7), ylim=c(-20,20),
#         qq.xlim=c(-7,5), col=2:3, do.par=FALSE,
#         qq.ylim=c(-100,100),
#         main="w1.ylim=c(-7,7) ylim=c(-20,20)\nqq.xlim=c(-7,5) qq.ylim=c(-100,100)")

par(org.par)
par(mfrow=c(3,3), mar=c(3, 3, 4, 1), mgp=c(2, 0.6, 0))

plotres(mod.glmnet.xmat, which=c(1,3,4), do.par=FALSE, #  w1.main="which=c(1,3,4)",
        w1.xlim=c(6,-6),
        w1.ylim=c(-5,5),
        w1.col=2:3,
        w1.main="TEST INTERACTION OF W1 ARGS PAGE 2\n\nwhich=c(1,3,4) w1.xlim=c(6,-6)\nw1.ylim=c(-5,5)) w1.col=2:3")

plotres(mod.glmnet.xmat, which=c(1,3,4), w1.cex.main=1,  do.par=FALSE, # w1.main="which=c(1,3,4)",
        xlim=c(-20,70),
        ylim=c(-60,60),
        w1.col=2:3,
        col=3:4,
        w1.main="which=c(1,3,4) ylim=c(-60,60))\nw1.col=2:3, col=3:4")

plotres(mod.glmnet.xmat, which=c(1,3,4), w1.cex.main=1, do.par=FALSE, # w1.main="which=c(1,3,4)",
        xlim=c(-20,70), w1.xlim=c(6,-6),
        ylim=c(-60,60), w1.ylim=c(-5,5),
        col=3:4,
        w1.main="which=c(1,3,4) xlim=c(9,-9), w1.xlim=c(6,-6)\nylim=c(-60,60), w1.ylim=c(-5,5)) w1.col=1:2, col=3:4")

par(org.par)

#-- make sure that we can work with all families

set.seed(2016)
par(mfrow=c(3,3), mar=c(3,3,3,1))
n <- 100
p <- 4
x <- matrix(rnorm(n*p), n, p)
g2 <- sample(1:2, n, replace=TRUE)
for(family in c("gaussian","binomial","poisson")) {
    mod <- glmnet(x,g2,family=family)
    plot(mod, xvar="lambda")
    plotres(mod, w1.xvar="lambda", main=paste("family", family),
            which=c(1,3), do.par=FALSE)
}
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
y <- cbind(time=ty, status=1-tcens) # y=Surv(ty,1-tcens) with library(survival)
glmnet.cox <- glmnet(x=x7, y=y, family="cox")
plot(glmnet.cox)
title("glmnet.cox", line=2)
plot_glmnet(glmnet.cox, xvar="norm")
plotres(glmnet.cox, which=3, do.par=FALSE)
par(org.par)

# test col argument
par(mfrow=c(2,3), mar=c(3,3,5,1), cex=1)
mod <- glmnet(as.matrix(mtcars[-1]), mtcars[,1])
plot_glmnet(mod, main="plot_glmnet  default")
plot_glmnet(mod, col=c(1,2,3,0,0,NA,0,0,0,0), main="col=c(1,2,3,0,0,NA,0,0,0,0)")
g <- "gray"
plot_glmnet(mod, col=c("black","red","green",g,g,g,g,g,"steelblue","darkorange"), main="col=c('black','red','green',g,g,g,g,g,'steelblue','darkorange')")
plot_glmnet(mod, col=c("black","red","green",0,0,0,0,0,"steelblue","darkorange"), main="col=c('black','red','green',0,0,0,0,0,'steelblue','darkorange')")
plot_glmnet(mod, col=c("black","red", 0), main="col=c('black','red', 0)") # test recycling, including 0
par(org.par)

source("test.epilog.R")
