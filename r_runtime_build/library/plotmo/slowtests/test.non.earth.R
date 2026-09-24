# test.non.earth.R: test plotmo on non-earth models
# Stephen Milborrow, Basley KwaZulu-Natal Mar 2011

source("test.prolog.R")
library(plotmo)
library(earth)
data(ozone1)
data(etitanic)
dopar <- function(nrows, ncols, caption = "")
{
    cat("                             ", caption, "\n")
    par(mfrow=c(nrows, ncols))
    par(oma = c(0, 0, 3, 0))
    par(mar = c(3, 3, 1.7, 0.5))
    par(mgp = c(1.6, 0.6, 0))
    par(cex = 0.7)
}
caption <- "test lm(log(doy) ~ vh+wind+humidity+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
set.seed(2020)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3, pt.pch=20, smooth.col="indianred",
       trace=2)
termplot(a)
par(org.par)

caption <- "test lm(log(doy) ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,5,caption)
a <- lm(log(doy) ~ vh + wind + humidity + temp + log(ibh), data=ozone1)
set.seed(2020)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.resp=3, pt.pch=20, clip=FALSE, smooth.col="indianred")
termplot(a)
par(org.par)

caption <- "test lm(doy ~ (vh+wind+humidity)^2, data=ozone1)"
dopar(4,3,caption)
a <- lm(doy ~ (vh+wind+humidity)^2, data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL)
# termplot(a) # termplot fails with Error in `[.data.frame`(mf, , i): undefined columns selected
par(org.par)

caption <- "test lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)"
dopar(4,3,caption)
a <- lm(doy^2 ~ vh+wind+humidity+I(wind*humidity)+temp+log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL)
termplot(a) # termplot draws a funky second wind plot
par(org.par)

caption <- "test lm with data=ozone versus attach(ozone)"
dopar(4,3,caption)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, degree1=c(1,2,4,5))
attach(ozone1)
a <- lm(log(doy) ~ I(vh*wind) + wind + I(humidity*temp) + log(ibh))
plotmo(a, do.par=FALSE, degree1=c(1,2,4,5))
detach(ozone1)
par(org.par)

# commented out because "$" in names is not yet supported
# a <- lm(log(ozone1$doy) ~ I(ozone1$vh*ozone1$wind) + log(ozone1$ibh))
# plotmo(a)

set.seed(1)
caption <- "test lm and glm a900..a902: damage~temp family=binomial data=orings"
dopar(2,3,caption)
library(faraway)
data(orings)
a900 <- lm(I(damage/6) ~ temp, data=orings)
plotmo(a900, do.par=FALSE, caption=caption, col.response=2, nrug=-1,
    main="lm(damage/6~temp)", smooth.col="indianred", trace=0)
response <- cbind(orings$damage, 6-orings$damage)
a901 <- glm(response ~ temp, family="binomial", data=orings)
set.seed(2020)
plotmo(a901, do.par=FALSE, col.response=2, nrug=-1,
       main="glm(response~temp)", smooth.col="indianred", trace=2)
a902 <- glm(cbind(damage, 6-damage)~temp, family="binomial", data=orings)
set.seed(2020)
plotmo(a902, do.par=FALSE, col.response=2, nrug=TRUE,
       main="glm(cbind(damage,6-damage)~temp)", trace=0)
termplot(a902, main="termplot")
plotmo(a902, type="link", main="type=\"link\"", do.par=F)
set.seed(2020)
plotmo(a902, type="response", main="type=\"response\"", col.response=2, do.par=F)
par(org.par)

set.seed(1)
caption <- "test glm(lot2~log(u),data=clotting,family=Gamma)"
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
set.seed(2020)
plotmo(a, do.par=FALSE, caption=caption, col.response=3, clip=FALSE, nrug=-1)
termplot(a)
plotmo(a, type="link", caption=paste("type=\"link\"", caption))
par(org.par)

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test plot.gam, with mgcv::gam(y ~ s(x) + s(x,z)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z1 <- runif(n);
y <- test1(cbind(x,z1)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z1))
set.seed(2020)
plotmo(a, do.par=FALSE, type2="contour", caption=caption,
      col.response=3, smooth.col="indianred",
      func=test1, func.col="indianred", func.lwd=5, func.lty=2, smooth.lwd=3)

plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", ylim=NA)
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)
n<-400
sig<-2
set.seed(1)
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f <- f0(x0) + f1(x1) + f2(x2)
e <- rnorm(n, 0, sig)
y <- f + e
test.func <- function(x) f0(x[,1]) + f1(x[,2]) + f2(x[,3])
library(mgcv)
caption <- "test mgcv::gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20)) (and extra persp plot)"
dopar(3,3,caption)
a <- gam(y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20))
plot(a, select=2)
plot(a, select=3)
plot(a, select=1)
plotmo(a, do.par=FALSE, type2="contour", caption=caption, xlab=NULL, main="", func=test.func, ngrid2=10, contour.drawlabels=FALSE)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, persp.the=-35)
par(org.par)

set.seed(1)
caption <- "test plot.gam, with mgcv::gam(doy~s(wind)+s(humidity,wind)+s(vh)+temp,data=ozone1)"
dopar(3,3,caption)
a <- gam(doy ~ s(wind) + s(humidity,wind) + s(vh) + temp, data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, type2="contour", degree1=c("wind","vh"), swapxy=T, xlab=NULL, main="", clip=FALSE)
plot(a, select=1)
plot(a, select=3)
plot(a, select=2)
plot(a, select=4)
par(org.par)

detach("package:mgcv")
library(gam)
caption <- "test gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
set.seed(1)
dopar(3,2,caption)
data(airquality)
airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
set.seed(2020)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3)
# termplot gives fishy looking wind plot, plotmo looks ok
# termplot(a) #TODO this fails with R2.5: dim(data) <- dim: attempt to set an attribute on NULL
detach("package:gam")
par(org.par)

library(mda)
caption <- "test mars and earth (expect not a close match)"
dopar(6,3,caption)
a <- mars( ozone1[, -1], ozone1[,1], degree=2)
b <- earth(ozone1[, -1], ozone1[,1], degree=2)
# this also tests trace=2 on a non formula model
plotmo(a, do.par=FALSE, caption=caption, trace=2)
plotmo(b, do.par=FALSE)
par(org.par)

caption <- "test mars and mars.to.earth(mars) (expect no degree2 for mars)"
dopar(6,3,caption)
a <- mars(ozone1[, -1], ozone1[,1], degree=2)
b <- mars.to.earth(a)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA)
plotmo(b, do.par=FALSE, ylim=NA)
par(org.par)

# check fix for bug reported by Martin Maechler:
# form <- Volume ~ .; a <- earth(form, data = trees); plotmo(a) fails

dopar(4,4, "test f <- O3 ~ .; a <- earth(f, data=ozone1)")
fa <- log(O3) ~ .
a <- earth(fa, data=ozone1, degree=2)
print(summary(a))
plot(a, do.par=FALSE)
set.seed(2020)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=c(1,2), col.response = "pink", smooth.col="indianred")
a <- lm(log(doy) ~ I(vh*wind) + I(humidity*temp) + log(ibh), data=ozone1)
plotmo(a, do.par=FALSE, degree1=1:2)
fa <- log(doy) ~ I(vh*wind) + I(humidity*temp) + log(ibh)
a <- lm(fa, data=ozone1)
plotmo(a, do.par=FALSE, degree1=1:2)
par(org.par)

# test inverse.func and func

caption <- "test inverse.func=exp"
a <- lm(log(Volume) ~ Girth + Height + I(Girth*Height), data=trees)
my.func <- function(x) -60 + 5 * x[,1] + x[,2] / 3
set.seed(2020)
plotmo(a, caption=caption, inverse.func = exp, col.response = "pink", func=my.func, func.col="gray", ngrid1=1000, type2="p", smooth.col="indianred")
par(org.par)

# se testing

caption = "level=.95, lm(doy~., data=ozone1) versus termplot"
dopar(6,3,caption)
a <- lm(doy~., data=ozone1)
plotmo(a, level=.95, do.par=FALSE, caption=caption)
termplot(a, se=2)
par(org.par)

caption <- "test different se options, level=.95, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(4,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95, level.shade="pink", level.shade2=3)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95, level.shade=3)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL, level=.95, level.shade=3)
par(org.par)

caption <- "test level=.95, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(2,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
termplot(a, se=2)
par(org.par)

caption <- "test level=.95 and inverse.func, lm(log(doy)~vh+wind+log(humidity),data=ozone1)"
dopar(3,3,caption)
a <- lm(log(doy) ~ vh + wind + log(humidity), data=ozone1)
plotmo(a, do.par=FALSE, caption=caption, ylim=NA, level=.95)
plotmo(a, do.par=FALSE, caption=caption, ylim=NULL, level=.95, inverse.func=exp)
termplot(a, se=2)
par(org.par)

caption <- "test level=.95, glm(lot2~log(u),data=clotting,family=Gamma)"
set.seed(1)
dopar(2,2,caption)
u = c(5,10,15,20,30,40,60,80,100)
lota = c(118,58,42,35,27,25,21,19,18)
clotting <- data.frame(u = u, lota = lota)
a <- glm(lota ~ log(u), data=clotting, family=Gamma)
set.seed(2020)
plotmo(a, do.par=FALSE, caption=caption, col.response=4, pt.pch=7, clip=FALSE, nrug=-1, level=.95, smooth.col="indianred")
termplot(a, se=2)
par(org.par)

if(length(grep("package:gam", search())))
    detach("package:gam")
library(mgcv)
set.seed(1)
caption <- "test level=.95, plot.gam, with mgcv::gam(y ~ s(x) + s(x,z1)) with response and func (and extra image plot)"
dopar(3,2,caption)
par(mar = c(3, 5, 1.7, 0.5))    # more space for left and bottom axis
test1 <- function(x,sx=0.3,sz=0.4)
    (pi**sx*sz)*(1.2*exp(-(x[,1]-0.2)^2/sx^2-(x[,2]-0.3)^2/sz^2)+
    0.8*exp(-(x[,1]-0.7)^2/sx^2-(x[,2]-0.8)^2/sz^2))
n <- 100
set.seed(1)
x <- runif(n);
z1 <- runif(n);
y <- test1(cbind(x,z1)) + rnorm(n) * 0.1
a <- gam(y ~ s(x) + s(x,z1))
set.seed(2020)
plotmo(a, do.par=FALSE, type2="contour", caption=caption, col.response=3, func=test1, func.col="magenta", level=.95)
plotmo(a, do.par=FALSE, degree1=F, degree2=1, type2="image", image.col=topo.colors(10),
        ylim=NA, level=.95, main="topo.colors")
plot(a, select=1)
plot(a, select=2)
plot(a, select=3)
par(org.par)

# TODO Following commented out because it causes:
#      Error: gam objects in the "gam" package do not support confidence intervals on new data
# detach("package:mgcv")
# library(gam)
# set.seed(1)
# caption <- "test level=.95, gam:gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind, Temp),data=airquality)"
# dopar(3,2,caption)
# data(airquality)
# airquality <- na.omit(airquality)   # plotmo doesn't know how to deal with NAs yet
# a <- gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data = airquality)
# set.seed(2020)
# plotmo(a, do.par=FALSE, caption=caption, ylim=NA, col.response=3, level=.95)
# # termplot(a)  #TODO this fails with R2.5: dim(data) <- dim: attempt to set an attribute on NULL
# detach("package:gam")
# par(org.par)

# test factors by changing wind to a factor

ozone2 <- ozone1
ozone2[,"wind"] <- factor(ozone2[,"wind"], labels=c(
    "wind0", "wind2", "wind3", "wind4", "wind5", "wind6",
    "wind7", "wind8", "wind9", "wind10", "wind11"))

# commented out because factors are not yet supported by plotmo.earth
# caption <- "test wind=factor, earth(O3 ~ ., data=ozone2)"
# a <- earth(doy ~ ., data=ozone2)
# set.seed(1)
# dopar(4,3,caption)
# set.seed(2020)
# plotmo(a, col.response="gray", level=.95, nrug=-1, do.par=FALSE, caption=caption)
# termplot(a)
# par(org.par)

caption <- "test wind=factor, lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone2)"
a <- lm(doy ~ vh + wind + I(humidity*temp) + log(ibh), data=ozone2)
set.seed(1)
dopar(4,3,caption)
plotmo(a, col.response="gray", level=.95, nrug=-1, do.par=FALSE, caption=caption, smooth.col="indianred")
termplot(a, se=2)
par(org.par)

caption <- "test level options"
dopar(2,2,caption)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade=0, caption=caption)
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade="orange")
plotmo(a, do.par=FALSE, degree1=2, degree2=FALSE, level=.95, level.shade2=0)
par(org.par)

caption <- "test wind=factor, glm(y ~ i + j, family=poisson())"
y <- c(18,17,15,20,10,20,25,13,12)
i <- gl(3,1,9)
j <- gl(3,3)
a <- glm(y ~ i + j, family=poisson())
set.seed(1)
dopar(2,2,caption)
plotmo(a, do.par=F, level=.95, nrug=1, caption=caption)
termplot(a, se=1, rug=T)
par(org.par)

if(length(grep("package:gam", search())))
   detach("package:gam")
caption <- "test wind=factor, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)"
library(mgcv)
a <- gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)
plotmo(a, level=.95, caption=caption)
caption <- "test wind=factor, clip=TRUE, gam(doy ~ vh + wind + s(humidity) + s(vh) + temp, data=ozone2)"
plotmo(a, level=.95, caption=caption, clip=FALSE)
# termplot doesn't work here so code commented out
# dopar(3,3,caption)
# plotmo(a, do.par=FALSE)
# termplot(a)
par(org.par)

# test lda and qda, and also col.response, pt.pch, and jitter
library(MASS)
etitanic2 <- etitanic
etitanic2$pclass <- as.numeric(etitanic$pclass)
etitanic2$sex <- as.numeric(etitanic$sex)
etitanic2$sibsp <- NULL
etitanic2$parch <- NULL
lda.model <- lda(survived ~ ., data=etitanic2)
set.seed(7)
plotmo(lda.model, caption="lda", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="posterior", nresponse=1, smooth.col="indianred",
       all2=TRUE, type2="image")
set.seed(8)
plotmo(lda.model, caption="lda with no jitter", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="posterior", nresponse=1,
       all2=TRUE, type2="image", jitter=0)
qda.model <- qda(survived ~ ., data=etitanic2)
set.seed(9)
plotmo(qda.model, caption="qda", clip=F,
       col.response=as.numeric(etitanic2$survived)+2, type="post", nresponse=2, smooth.col="indianred",
       all2=TRUE, type2="image", jitter.resp=.6, pch.resp=20)

# test plotmo.y from the 2nd argument of the model function (non-formula interface)
lcush <- data.frame(Type=as.numeric(Cushings$Type), log(Cushings[,1:2]))[1:21,]
a <- qda(lcush[,2:3], lcush[,1])
set.seed(2020)
plotmo(a, type="class", all2=TRUE,
       caption= "plotmo.y from 2nd argument of call (qda)",
       type2="contour", ngrid2=100, contour.nlevels=2, contour.drawlabels=FALSE,
       col.response=as.numeric(lcush$Type)+1,
       pt.pch=as.character(lcush$Type))
par(org.par)

# # example from MASS (works, but removed because unnecessary test)
# predplot <- function(object, main="", len = 100, ...)
# {
#     plot(Cushings[,1], Cushings[,2], log="xy", type="n",
#          xlab = "Tetrahydrocortisone", ylab = "Pregnanetriol", main = main)
#     for(il in 1:4) {
#         set <- Cushings$Type==levels(Cushings$Type)[il]
#         text(Cushings[set, 1], Cushings[set, 2],
#              labels=as.character(Cushings$Type[set]), col = 2 + il) }
#     xp <- seq(0.6, 4.0, length=len)
#     yp <- seq(-3.25, 2.45, length=len)
#     cushT <- expand.grid(Tetrahydrocortisone = xp,
#                          Pregnanetriol = yp)
#     Z <- predict(object, cushT, ...); zp <- as.numeric(Z$class)
#     zp <- Z$post[,3] - pmax(Z$post[,2], Z$post[,1])
#     contour(exp(xp), exp(yp), matrix(zp, len),
#             add = TRUE, levels = 0, labex = 0)
#     zp <- Z$post[,1] - pmax(Z$post[,2], Z$post[,3])
#     contour(exp(xp), exp(yp), matrix(zp, len),
#             add = TRUE, levels = 0, labex = 0)
#     invisible()
# }
# par(mfrow=c(2,2))
# cush <- log(as.matrix(Cushings[, -3]))
# tp <- Cushings$Type[1:21, drop = TRUE]
# set.seed(203)
# cush.data <- data.frame(tp, cush[1:21,])
# a <- qda(tp~., data=cush.data)
# predplot(a, "QDA example from MASS")
# set.seed(2020)
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1)
# set.seed(2020)
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1, drawlabels=F, nlevels=2)
# set.seed(2020)
# plotmo(a, type="class", all2=TRUE, type2="contour", degree1=NA, do.par=FALSE,
#        col.response=as.numeric(cush.data$tp)+1, drawlabels=F, nlevels=2, ngrid2=100)
# par(org.par)

library(rpart)
data(kyphosis)
# kyphosis data, earth model
a <- earth(Kyphosis ~ ., data=kyphosis, degree=2, glm=list(family=binomial))
cat("summary(a): (Kyphosis)\n")
print(summary(a))
par(mfrow=c(3, 3))
par(mar=c(3, 3, 2, .5))  # small margins to pack figs in
set.seed(9) # for jitter
set.seed(2020)
plotmo(a, do.par=F, type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       clip=F)
plotmo(a, do.par=F, clip=F, degree1=0)
par(org.par)

# kyphosis data, rpart models (also test ngrid2)
fit1 <- rpart(Kyphosis ~ ., data=kyphosis)
plotres(fit1, SHOWCALL=TRUE)
par(mfrow=c(3, 3))
par(mar=c(.5, 0.5, 2, .5), mgp = c(1.6, 0.6, 0))  # b l t r small margins to pack figs in
library(rpart.plot)
prp(fit1, main="rpart kyphosis\nno prior")
plotmo(fit1, degree1=NA, do.par=F, main="", persp.theta=220, nresponse=2)
par(mar=c(4, 4, 2, .5))
set.seed(2020)
plotmo(fit1, nresp=2, degree1=FALSE, do.par=F, main="", type2="image", # test default type="prob"
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=ifelse(kyphosis$Kyphosis=="present", "p", "a"),
       image.col=gray(10:4/10), ngrid2=30)
par(mar=c(.5, 0.5, 2, .5))  # b l t r small margins to pack figs in
plotmo(fit1, type="class", degree1=NA, do.par=F, main="type=\"class\"")
# with type="prob" and response has two columns,
# nresponse should automatically default to column 2
plotmo(fit1, type="prob", degree1=0, do.par=F, main="type=\"prob\"",
       clip=F, ngrid2=50, persp.border=NA, trace=1)
set.seed(2020)
plotmo(fit1, type="prob", nresp=2, degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10), ngrid2=5)
# better rpart model with prior
fit2 <- rpart(Kyphosis ~ ., data=kyphosis, parms=list(prior=c(.65,.35)))
prp(fit2, main="rpart kyphosis\nwith prior, better model")
plotmo(fit2, type="v", degree1=NA, do.par=F, main="", persp.theta=220, ngrid2=10)
par(mar=c(4, 4, 2, .5))
set.seed(2020)
plotmo(fit2, type="v", degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10), ngrid2=100)
par(org.par)

plotmo(fit1, type="prob", nresponse=1, persp.border=NA, persp.col="pink", all1=TRUE, all2=TRUE,
       caption="plotmo rpart fit1, all1=TRUE, all2=TRUE")
expect.err(try(plotmo(fit1, type="none.such1")))

# rpart model with ozone data
data(ozone1)
par(mfrow=c(4,4))
par(mar=c(.5, 0.5, 2, .5), cex=.6, mgp = c(1.6, 0.6, 0))  # b l t r small margins to pack figs in
a1 <- rpart(O3~temp+humidity, data=ozone1)
prp(a1, main="rpart model with ozone data\n(temp and humidity only)\n")
plotmo(a1, do.par=F, degree1=0, main="rpart", persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a1, type="class")))
# compare to a linear and earth model
a3 <- lm(O3~temp+humidity, data=ozone1)
plotmo(a3, do.par=F, clip=F, main="lm", degree1=0, all2=TRUE, persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a3, type="none.such2")))
a <- earth(O3~temp+humidity, data=ozone1, degree=2)
plotmo(a, do.par=F, clip=F, main="earth", degree1=NA, persp.ticktype="detail", persp.nticks=2)
expect.err(try(plotmo(a, type="none.such3")))
expect.err(try(plotmo(a, type=c("abc", "def"))))
par(org.par)

# detailed rpart model
par(mfrow=c(3,3))
a1 <- rpart(O3~., data=ozone1)
prp(a1, cex=.9, main="rpart model with full ozone data")
plotmo(a1, type="vector", do.par=F, degree1=NA, persp.ticktype="detail",
       persp.nticks=3, degree2=2:3)
par(org.par)

plotmo(a1, persp.border=NA, all1=TRUE, all2=TRUE,
       caption="plotmo rpart a1, all1=TRUE, all2=TRUE")

library(tree)
tree1 <- tree(O3~., data=ozone1)
plotmo(tree1)
plotres(tree1)

# rpart data with NAs

rpart.airquality <- rpart(Ozone~., data=airquality) # airquality has NAs in response and variables
plotmo <- plotmo(rpart.airquality, trace=0, SHOWCALL=TRUE)
print(rpart.rules(rpart.airquality))

airquality.nonaOzone <- subset(airquality, !is.na(Ozone)) # no NAs in response but NAs in variables
rpart.nonaOzone <- rpart(Ozone~., data=airquality.nonaOzone)
print(rpart.rules(rpart.nonaOzone))
plotmo.nonaOzone <- plotmo(rpart.nonaOzone, trace=0, SHOWCALL=TRUE)
airquality.nonaOzone$Ozone <- NULL
stopifnot(identical(plotmo.nonaOzone, airquality.nonaOzone))

# test xflip and yflip

par(mfrow=c(4, 4))
par(mgp = c(1.6, 0.6, 0))
par(mar=c(4, 4, 2, .5))

flip.test1 <- rpart(Kyphosis ~ ., data=kyphosis)
set.seed(2020)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10))
set.seed(2020)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="xflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       xflip=T)
set.seed(2020)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="yflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       yflip=T)
set.seed(2020)
plotmo(flip.test1, type="prob", nresp=2, degree1=NA, do.par=F, main="xflip and yflip", type2="image",
       col.response=ifelse(kyphosis$Kyphosis=="present", "red", "lightblue"),
       pt.pch=20, image.col=gray(10:4/10),
       xflip=T, yflip=T)

flip.test2 <- earth(O3~., data=ozone1, degree=2)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="", type2="cont")
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip", type2="cont",
       xflip=T)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="yflip", type2="cont",
       yflip=T)
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip and yflip", type2="cont",
       xflip=T, yflip=T)

cat("Expect warnings: ignoring xflip=TRUE for persp plot\n")
plotmo(flip.test2, degree1=NA, degree2=2, do.par=F, main="xflip and yflip", type2="persp",
       xflip=T, yflip=T)

library(randomForest)
data(etitanic)
etit <- etitanic[1:300,]

cat("=== rf.regression ===\n")
set.seed(2016)
# Expect Warning: The response has five or fewer unique values.  Are you sure you want to do regression?
rf.regression <- randomForest(survived~., data=etit, ntree=100, importance = FALSE)
plotmo(rf.regression, trace=1)

cat("=== rf.regression.importance ===\n")
set.seed(2016)
# Expect Warning: The response has five or fewer unique values.  Are you sure you want to do regression?
rf.regression.importance <- randomForest(survived~., data=etit, ntree=100, importance = TRUE)
plotmo(rf.regression.importance, trace=1)

etit <- etitanic[1:300,]
etit$survived <- factor(ifelse(etit$survived == 1, "survived", "died"),
                       levels = c("survived", "died"))
cat("=== rf.classification ===\n")
set.seed(2016)
rf.classification <- randomForest(survived~., data=etit, ntree=100, importance = FALSE)
plotmo(rf.classification, trace=1, type="prob", nresponse="surv", SHOWCALL=TRUE)
plotmo(rf.classification, trace=1, type="prob", nresponse="died", degree2=0, SHOWCALL=TRUE)

cat("=== rf.classification.importance ===\n")
set.seed(2016)
rf.classification.importance <- randomForest(survived~., data=etit, ntree=100, importance = TRUE)
plotmo(rf.classification.importance, trace=1, type="prob", nresponse="surv", SHOWCALL=TRUE)

cat("=== plotres randomForest ===\n")
plotres(rf.regression)
plotres(rf.regression.importance)
# TODO residuals are in range 0 to 1
plotres(rf.classification, type="prob", nresponse="surv")
plotres(rf.classification.importance, type="prob", nresponse="surv")

#--- fda ------------------------------------------------------------------------------

par(org.par)
par(mfrow=c(4,5))
par(mar = c(3, 2, 3, .1)) # b, l, t, r
par(mgp = c(1.5, .5, 0))
fda.earth <- fda(Species~., data=iris, keep.fitted=TRUE, method=earth, keepxy=TRUE)
fda.polyreg <- fda(Species~., data=iris, keep.fitted=TRUE, keepxy=TRUE)
fda.bruto <- fda(Species~., data=iris, keep.fitted=TRUE, method=bruto)

# 'fda.polyreg$fit' does not have a 'call' field or 'x' and 'y' fields
expect.err(try(plotmo(fda.polyreg$fit, type="variates", nresponse=1, clip=F, do.par=F)))

plot(1, main="plotmo with fda", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.2, xpd=NA)

plotmo(fda.earth, type="variates", nresponse=1, clip=F, do.par=F)

plot(1, main="plotmo with fda.earth$fit", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.2, xpd=NA)

plotmo(fda.earth$fit, nresponse=1, clip=F, do.par=F)

plot(1, main="", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.5, xpd=NA)

plot(fda.earth)
plotmo(fda.earth, clip=F, do.par=F) # default type is class

plot(fda.polyreg)
plotmo(fda.polyreg, type="variates", nresponse=1, clip=F, do.par=F, degree1=c(1,3,4))
plot(1, main="", xaxt="n", yaxt="n", xlab="", ylab="",
     type="n", bty="n", cex.main=1.5, xpd=NA)

par(mfrow=c(3,3))
par(mar = c(3, 2, 3, .1)) # b, l, t, r
par(mgp = c(1.5, .5, 0))
plot(fda.bruto)
plotmo(fda.bruto, type="variates", nresponse=1, do.par=F)
par(org.par)

# neural net package
# for speed we use artificial data because neuralnet is very slow on say trees
library(neuralnet)
n <- 20
set.seed(3)
x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1) # x2 is noise
y <- x1^2
data <- data.frame(y=y, x1=x1, x2=x2)
colnames(data) <- c("y","x1", "x2")
set.seed(3)
nn <- neuralnet(y~x1+x2, data=data, hidden=3, rep=3)
print(head(plotmo:::predict.nn(nn, rep="best", trace=TRUE)))
set.seed(2020)
plotmo(nn, trace=1, col.response=2, all2=TRUE, SHOWCALL=TRUE)
# trace=0 below to test hushing of message "assuming "y" in the model.frame is the response, because object$terms is NULL"
set.seed(2020)
plotmo(nn, trace=0, col.response=2, predict.rep="best", SHOWCALL=TRUE)
plotres(nn, trace=0, info=TRUE, SHOWCALL=TRUE)
plotres(nn, trace=1, info=TRUE, predict.rep="best", SHOWCALL=TRUE)

library(nnet)
data(iris3)
set.seed(301)
samp <- c(sample(1:50,25), sample(51:100,25), sample(101:150,25))
ird <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                  species=factor(c(rep("seto",50), rep("vers", 50), rep("virg", 50))))
ir.nn2 <- nnet(species ~ ., data = ird, subset = samp, size = 2, rang = 0.1,
               decay = 5e-4, maxit = 20)
plotmo(ir.nn2, nresponse=1, type="class", all2=T, degree2=2:6)
plotmo(ir.nn2, nresponse=2, clip=F, all2=T, degree2=1:5)
plotres(ir.nn2, nresponse=2)

library(biglm)
data(trees)
ff <- log(Volume)~log(Girth)+log(Height)
chunk1 <- trees[1:20,]
chunk2 <- trees[20:31,]
biglm <- biglm(ff,chunk1)
biglm <- update(biglm, chunk2)
plotmo(biglm, pt.col=2, SHOWCALL=TRUE)
plotres(biglm, SHOWCALL=TRUE)

library(adabag)
data(iris)
set.seed(2015)
# mfinal=3 for speed during testing
mod.boosting <- boosting(Species~., data=iris, mfinal=3)
mod.bagging <- bagging(Species~., data=iris, mfinal=3)
dopar(4, 4, caption="adabag package")
plotmo(mod.boosting, nresponse=1, ylim=c(0,1), do.par=FALSE) # default type="prob"
plotmo(mod.boosting, type="class", do.par=FALSE)
plotmo(mod.bagging, nresponse=1, ylim=c(0,1), do.par=FALSE)
plotmo(mod.bagging, nresponse=1, type="votes", do.par=FALSE)
par(org.par)

library(e1071)
data(iris)
x.iris <- subset(iris, select=-Species)
y.iris <- iris$Species
set.seed(2016)
svm.xy <- svm(x.iris, y.iris, probability=FALSE)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
expect.err(try(plotmo(svm.xy, prob=TRUE, nresponse="vers", do.par=TRUE, all2=TRUE))) # probability=FALSE in call to svm
plotmo(svm.xy, decision=TRUE,
       nresponse="vers", do.par=FALSE, all2=TRUE, degree1=2:3, degree2=6)
svm.xy <- svm(x.iris, y.iris, probability=TRUE)
plotmo(svm.xy, prob=TRUE,
       nresponse="vers", do.par=FALSE, all2=TRUE, degree1=2:3, degree2=6)
set.seed(2016)
svm.form <- svm(Species ~ ., data=iris, probability=T)
plotmo(svm.form, predict.p=TRUE,
       nresponse="vers", do.par=FALSE, all2=TRUE, degree1=2:3, degree2=6)
expect.err(try(plotmo(svm.form, decision.values=TRUE, probab=TRUE))) # not both
plotres(svm.form, predict.prob=TRUE, nresponse="vers", info=TRUE)
plotres(svm.form, jitter=5, info=TRUE)
par(org.par)

source("test.epilog.R")
