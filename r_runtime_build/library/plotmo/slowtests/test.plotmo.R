# test.plotmo.R: regression tests for plotmo
# Stephen Milborrow, Petaluma Jan 2007

print(R.version.string)

source("test.prolog.R")
library(earth)
options(warn=1) # print warnings as they occur
data(etitanic)
make.space.for.caption <- function(caption="CAPTION")
{
    oma <- par("oma")
    needed <- 3
    # adjust for newlines in caption
    newlines <- grep("\n", caption)
    if(length(newlines) > 0)
        needed <- needed + .5 * newlines # .5 seems enough although 1 in theory
    if(!is.null(caption) && any(nchar(caption)) && oma[3] <= needed) {
        oma[3] <- needed
        par(oma=oma)
    }
}
dopar <- function(nrows, ncols, caption = "")
{
    cat("                             ", caption, "\n")
    make.space.for.caption(caption)
    par(mfrow=c(nrows, ncols))
    par(mar = c(3, 3, 1.7, 0.5))
    par(mgp = c(1.6, 0.6, 0))
    par(cex = 0.7)
}
example(plotmo)
caption <- "basic earth test of plotmo"
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, degree1=2, degree2=4, caption=caption, trace=-1)

caption <- "test 5 x 5 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=51, pmethod="n", degree=2)
plotmo(a, caption=caption, trace=1)

caption <- "test 4 x 4 layout with ylab"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=30, pmethod="n", degree=2)
plotmo(a, caption=caption, trace=2)

caption <- "test 3 x 3 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=16, pmethod="n", degree=2)
plotmo(a, caption=caption, trace=3)

caption <- "test 2 x 2 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=9, pmethod="n", degree=2)
plotmo(a, caption=caption)

caption <- "test 1 x 1 layout"
dopar(1,1,caption)
a <- earth(O3 ~ ., data=ozone1, nk=4, pmethod="n", degree=2)
plotmo(a, caption=caption)

caption <- "test plotmo basic params"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(3,2,caption)
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, caption=caption,
        main="test main", xlab="test xlab", ylab="test ylab")
plotmo(a, do.par=FALSE, degree1=F, degree2=4, grid.func=mean, persp.col="white", ngrid2=10, persp.phi=40)
set.seed(2016)
plotmo(a, do.par=FALSE, degree1=1, degree1.lty=2, degree1.lwd=4, degree1.col=2, nrug=TRUE, degree2=F, main="nrug=300")
plotmo(a, do.par=FALSE, degree1=1, nrug=-1, degree2=F, main="nrug=TRUE")
set.seed(2016)
plotmo(a, do.par=FALSE, degree1=1, nrug=10, ngrid1=50, degree2=F, main="ngrid1=50 nrug=10")
plotmo(a, do.par=FALSE, degree1=NA, degree2=1, persp.phi=60) # graph args

caption <- "test plotmo xlim and ylim"
a <- earth(O3 ~ ., data=ozone1, degree=2)
dopar(5,3,caption)
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, caption=caption, xlab="ylim=default")
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=NA, xlab="ylim=NA")
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, ylim=c(0,20), xlab="ylim=c(0,20)")
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, xlim=c(190,250), xlab="xlim=c(190,250)")
plotmo(a, do.par=FALSE, degree1=2:3, degree2=4, xlim=c(190,250), ylim=c(11,18), xlab="xlim=c(190,250), ylim=c(11,18)")

# check various types of predictors with grid.func and ndiscrete

varied.type.data <- data.frame(
    y    = 1:13,
    num  = c(1, 3, 2, 3, 4, 5, 6, 4, 5, 6.5, 3, 6, 5), # 7 unique values (but one is non integral)
    int  = c(1L, 1L, 3L, 3L, 4L, 4L, 3L, 5L, 3L, 6L, 7L, 8L, 10L), # 8 unique values
    bool = c(F, F, F, F, F, T, T, T, T, T, T, T, T),
    date = as.Date(
           c("2018-08-01", "2018-08-02", "2018-08-03",
             "2018-08-04", "2018-08-05", "2018-08-06",
             "2018-08-07", "2018-08-08", "2018-08-08",
             "2018-08-08", "2018-08-10", "2018-08-11", "2018-08-11")),
    ord  = ordered(c("ord3", "ord3", "ord3",
                     "ord1", "ord2", "ord3",
                     "ord1", "ord2", "ord3",
                     "ord1", "ord1", "ord1", "ord1"),
                   levels=c("ord1", "ord3", "ord2")),
    fac  = as.factor(c("fac1", "fac1", "fac1",
                       "fac2", "fac2", "fac2",
                       "fac3", "fac3", "fac3",
                       "fac1", "fac2", "fac3", "fac3")),
    str  = c("str1", "str1", "str1", # will be treated like a factor
             "str2", "str2", "str2",
             "str3", "str3", "str3",
             "str3", "str3", "str3", "str3"))

varied.type.lm <- lm(y ~ ., data = varied.type.data)
print(summary(varied.type.lm))
set.seed(2018)
plotres(varied.type.lm, info=TRUE)
plotmo(varied.type.lm, pmethod="apartdep", all2=TRUE, ticktype="d", col.response="red", caption="varied.type.lm\npmethod=\"apartdep\" default grid func")
plotmo(varied.type.lm, all2=TRUE, ticktype="d", col.response="red", caption="varied.type.lm\ndefault grid func")
plotmo(varied.type.lm, all2=TRUE, ndiscre=1, caption="varied.type.lm\nndiscrete=1")
plotmo(varied.type.lm, all2=TRUE, ndiscr=2, caption="varied.type.lm\nndiscrete=2")
plotmo(varied.type.lm, all2=TRUE, ndis=100, caption="varied.type.lm\nndiscrete=100")
cat("grid.func=median:\n")
plotmo(varied.type.lm, all2=TRUE, grid.func=median, caption="varied.type.lm\ngrid.func=median")
cat("grid.func=quantile:\n")
plotmo(varied.type.lm, all2=TRUE, grid.func=function(x, ...) quantile(x, 0.5), caption="varied.type.lm\ngrid.func=function(x, ...) quantile(x, 0.5)")
cat("grid.func=mean:\n")
plotmo(varied.type.lm, all2=TRUE, grid.func=mean, caption="varied.type.lm\ngrid.func=mean")

varied.type.earth <- earth(y ~ ., data = varied.type.data, thresh=0, penalty=-1, trace=1)
print(summary(varied.type.earth))
set.seed(2018)
plotres(varied.type.earth, info=TRUE)
plotmo(varied.type.earth, all1=TRUE, all2=TRUE, persp.ticktype="d", col.response="red")

# term.plot calls predict.earth with an se parameter, even with termplot(se=FALSE)

caption <- "basic earth test against termplot"
dopar(4,4,caption)
make.space.for.caption("test caption1")
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, do.par=FALSE, ylim=NA, caption=caption, degree2=FALSE)
cat("Ignore warning: predict.earth ignored argument \"se.fit\"\n")
termplot(a)

caption <- "test change order of earth predictors and cex"
dopar(4,4,caption)
# minspan=1 to force two degree2 graphs for the test (wasn't necessary in old versions of earth)
a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2, minspan=1)
plotmo(a, do.par=FALSE, ylim=NA, caption=caption, degree2=c(1,2), cex=1.2)
termplot(a)

caption <- "test all1=TRUE"
a <- earth(doy ~ humidity + temp + wind, data=ozone1, degree=2)
plotmo(a, caption=caption, all1=TRUE, persp.ticktype="d", persp.nticks=2)
caption <- "test all2=TRUE"
print(summary(a))
plotmo(a, caption=caption, all2=TRUE)

oz <- ozone1[150:200,c("O3","temp","humidity","ibh")]
a.glob <- earth(O3~temp+humidity, data=oz, degree=2)
ad.glob <- earth(oz[,2:3], oz[,1], degree=2)
func1 <- function()
{
    caption <- "test environments and finding the correct data"
    dopar(4,4,caption)
    set.seed(2016)

    plotmo(a.glob, do.par=FALSE, main="a.glob oz",
          degree1=1, all2=1, degree2=1, type2="im",
          col.response=3, pt.pch=20, trace=2)
    mtext(caption, outer=TRUE, font=2, line=1.5, cex=1)
    plotmo(ad.glob, do.par=FALSE, main="ad.glob oz",
          degree1=1, all2=1, degree2=1, type2="im",
          col.response=3, pch.response=20, trace=2) # pch.response test backcompat

    a <- earth(O3~temp+humidity, data=oz, degree=2)
    plotmo(a, do.par=FALSE, main="a oz",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

    ad <- earth(oz[,2:3], oz[,1], degree=2)
    plotmo(ad, do.par=FALSE, main="ad oz",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

    oz.org <- oz
    oz10 <- 10 * oz # multiply by 10 so we can see by the axis labels if right data is being used
    oz <- oz10      # oz is now local to this function, but multiplied by 10
    a.oz10 <- earth(O3~temp+humidity, data=oz, degree=2)
    a.oz10.keep <- earth(O3~temp+humidity, data=oz, degree=2, keepxy=TRUE)
    plotmo(a.oz10, do.par=FALSE, main="a oz10",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

    ad.oz10 <- earth(oz[,2:3], oz[,1], degree=2)
    ad.oz10.keep <- earth(oz[,2:3], oz[,1], degree=2, keepxy=TRUE)
    plotmo(ad.oz10, do.par=FALSE, main="ad oz10",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

    func2 <- function() {
        a.func <- earth(O3 ~ temp +   humidity, data=oz10, degree=2)
        plotmo(a.func, do.par=FALSE, main="a.func oz10",
               degree1=1, all2=1, degree2=1, type2="im",
               col.response=3, pt.pch=20)

        ad.func <- earth(oz10[,2:3], oz10[,1], degree=2)
        plotmo(ad.func, do.par=FALSE, main="ad.func oz10",
               degree1=1, all2=1, degree2=1, type2="im",
               col.response=3, pt.pch=20)

        caption <- "test environments and finding the correct data, continued"
        dopar(4,4,caption)

        oz <- .1 * oz.org
        a.func <- earth(O3~temp+ humidity , data=oz, degree=2)
        plotmo(a.func, do.par=FALSE, main="a.func oz.1",
               degree1=1, all2=1, degree2=1, type2="im",
               col.response=3, pt.pch=20)

        ad.func <- earth(oz[,2:3], oz[,1], degree=2)
        plotmo(ad.func, do.par=FALSE, main="ad.func oz.1",
               degree1=1, all2=1, degree2=1, type2="im",
               col.response=3, pt.pch=20)

        plotmo(a.oz10.keep, do.par=FALSE, main="func1:a.oz10.keep",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

        plotmo(ad.oz10.keep, do.par=FALSE, main="func1:ad.oz10.keep",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20)

        cat("Expect error msg: formal argument \"do.par\" matched by multiple actual arguments\n")
        expect.err(try(plotmo(a.oz10, do.par=FALSE, main="func1:a.oz10",
           degree1=1, all2=1, degree2=1, type2="im",
           col.response=3, pt.pch=20, do.par=FALSE)))
    }
    func2()

    y  <- 3:11
    x1 <- c(1,3,2,4,5,6,6,6,6)
    x2 <- c(2,3,4,5,6,7,8,9,10)
    frame <- data.frame(y=y, x1=x1, x2=x2)
    foo <- function()
    {
        lm.18.out <- lm(y~x1+x2, model=FALSE)
        x1[2] <- 18
        y[3] <- 19
        frame <- data.frame(y=y, x1=x1, x2=x2)
        list(lm.18.out   = lm.18.out,
             lm.18       = lm(y~x1+x2),
             lm.18.keep  = lm(y~x1+x2, x=TRUE, y=TRUE),
             lm.18.frame = lm(y~x1+x2, data=frame))
    }
    temp <- foo()
        lm.18.out   <- temp$lm.18.out
        lm.18       <- temp$lm.18
        lm.18.keep  <- temp$lm.18.keep
        lm.18.frame <- temp$lm.18.frame

    # following should all use the x1 and y inside foo

    cat("==lm.18.out\n")
    plotmo(lm.18.out, main="lm.18.out",
           do.par=FALSE, degree1=1, clip=FALSE, ylim=c(0,20),
           col.response=2, pt.pch=20)

    cat("==lm.18\n")
    plotmo(lm.18, main="lm.18",
           do.par=FALSE, degree1=1, clip=FALSE, ylim=c(0,20),
           col.response=2, pt.pch=20)

    cat("==lm.18.keep\n")
    plotmo(lm.18.keep, main="lm.18.keep", trace=2,
           do.par=FALSE, degree1=1, clip=FALSE, ylim=c(0,20),
           col.response=2, pt.pch=20)

    cat("==lm.18.frame\n")
    plotmo(lm.18.frame, main="lm.18.frame",
           do.par=FALSE, degree1=1, clip=FALSE, ylim=c(0,20),
           col.response=2, pt.pch=20)
}
func1()

caption <- "test earth formula versus x,y model"
# dopar(4,4,caption)
# mtext(caption, outer=TRUE, font=2, line=1.5, cex=1)
a <- earth(O3 ~ ., data=ozone1, degree=2)
plotmo(a, caption="test earth formula versus xy model (formula)")
a <- earth(ozone1[, -1], ozone1[,1], degree=2)
plotmo(a, caption="test earth formula versus xy model (xy)")

# single predictor
caption <- "test earth(O3~wind, data=ozone1, degree=2), single predictor"
dopar(2,2,caption)
a <- earth(O3~wind, data=ozone1, degree=2)
plotmo(a)

caption = "se=2, earth(doy~humidity+temp+wind, data=ozone1) versus termplot (expect no se lines)"
dopar(3,3,caption)
mtext(caption, outer=TRUE, font=2, line=1.5, cex=1)
# minspan=1 to force two degree2 graphs for the test (wasn't necessary in old versions of earth)
a <- earth(doy~humidity + temp + wind, data=ozone1, degree=2, minspan=1)
cat("Ignore warning: predict.earth ignored argument \"se\"\n")
termplot(a)
plotmo(a, do.par=FALSE, ylim=NA, degree2=c(1:2), clip=FALSE, caption=caption)

# test fix to bug reported by Joe Retzer, FIXED Dec 7, 2007
N <- 650
set.seed(2007)
q_4    <- runif(N, -1, 1)
q_2102 <- runif(N, -1, 1)
q_2104 <- runif(N, -1, 1)
q_3105 <- runif(N, -1, 1)
q_3106 <- runif(N, -1, 1)
q_4104 <- runif(N, -1, 1)
q_6101 <- runif(N, -1, 1)
q_6103 <- runif(N, -1, 1)
q_7104 <- runif(N, -1, 1)
q_3109 <- runif(N, -1, 1)
q_4103 <- runif(N, -1, 1)
q_2111 <- runif(N, -1, 1)
q_3107 <- runif(N, -1, 1)
q_3101 <- runif(N, -1, 1)
q_3104 <- runif(N, -1, 1)
q_7107 <- runif(N, -1, 1)
depIndex <- sin(1.0 * q_4 + rnorm(650, sd=.8)) + sin(1.8 * q_2102 + rnorm(650, sd=.8)) + sin(1.3 * q_2104 + rnorm(650, sd=.8)) + sin(1.4 * q_3105 + rnorm(650, sd=.8)) +
            sin(1.5 * q_3106 + rnorm(650, sd=.8)) + sin(1.6 * q_4104 + rnorm(650, sd=.8)) + sin(1.8 * q_6101 + rnorm(650, sd=.8)) + sin(1.8 * q_6103 + rnorm(650, sd=.8)) +
            sin(1.9 * q_7104 + rnorm(650, sd=.8)) + sin(2.0 * q_3109 + rnorm(650, sd=.8))

regDatCWD <- as.data.frame(cbind(depIndex, q_4, q_2102, q_2104, q_3105, q_3106, q_4104, q_6101, q_6103, q_7104, q_3109, q_4103, q_2111, q_3107, q_3101, q_3104, q_7107))
cat("--plotmo(earthobj5)--\n")
earthobj5 <- earth(depIndex ~  q_4+q_2102+q_2104+q_3105+q_3106+q_4104+q_6101+q_6103+q_7104+q_3109+q_4103+q_2111+q_3107+q_3101+q_3104+q_7107, data=regDatCWD)
print(summary(earthobj5, digits = 2))
plotmo(earthobj5)

# long predictor names

a.rather.long.in.fact.very.long.name.q_4 <- q_4
a.rather.long.in.fact.very.long.name.q_2102 <- q_2102
a.rather.long.in.fact.very.long.name.q_2104 <- q_2104
a.rather.long.in.fact.very.long.name.q_3105 <- q_3105
a.rather.long.in.fact.very.long.name.q_3106 <- q_3106
a.rather.long.in.fact.very.long.name.q_4104 <- q_4104
a.rather.long.in.fact.very.long.name.q_6101 <- q_6101
a.rather.long.in.fact.very.long.name.q_6103 <- q_6103
a.rather.long.in.fact.very.long.name.q_7104 <- q_7104
a.rather.long.in.fact.very.long.name.q_3109 <- q_3109
a.rather.long.in.fact.very.long.name.q_4103 <- q_4103
a.rather.long.in.fact.very.long.name.q_2111 <- q_2111
a.rather.long.in.fact.very.long.name.q_3107 <- q_3107
a.rather.long.in.fact.very.long.name.q_3101 <- q_3101
a.rather.long.in.fact.very.long.name.q_3104 <- q_3104
a.rather.long.in.fact.very.long.name.q_7107 <- q_7107
a.rather.long.in.fact.very.long.name.for.the.response <- depIndex
a.rather.long.in.fact.very.long.name.for.the.dataframe <-
        as.data.frame(cbind(
                a.rather.long.in.fact.very.long.name.for.the.response,
                a.rather.long.in.fact.very.long.name.q_4,
                a.rather.long.in.fact.very.long.name.q_2102,
                a.rather.long.in.fact.very.long.name.q_2104,
                a.rather.long.in.fact.very.long.name.q_3105,
                a.rather.long.in.fact.very.long.name.q_3106,
                a.rather.long.in.fact.very.long.name.q_4104,
                a.rather.long.in.fact.very.long.name.q_6101,
                a.rather.long.in.fact.very.long.name.q_6103,
                a.rather.long.in.fact.very.long.name.q_7104,
                a.rather.long.in.fact.very.long.name.q_3109,
                a.rather.long.in.fact.very.long.name.q_4103,
                a.rather.long.in.fact.very.long.name.q_2111,
                a.rather.long.in.fact.very.long.name.q_3107,
                a.rather.long.in.fact.very.long.name.q_3101,
                a.rather.long.in.fact.very.long.name.q_3104,
                a.rather.long.in.fact.very.long.name.q_7107))

cat("--a.rather.long.in.fact.very.long.name.for.the...A--\n")
a.rather.long.in.fact.very.long.name.for.the.modelA <-
        earth(a.rather.long.in.fact.very.long.name.for.the.response ~
                a.rather.long.in.fact.very.long.name.q_4 +
                a.rather.long.in.fact.very.long.name.q_2102 +
                a.rather.long.in.fact.very.long.name.q_2104 +
                a.rather.long.in.fact.very.long.name.q_3105 +
                a.rather.long.in.fact.very.long.name.q_3106 +
                a.rather.long.in.fact.very.long.name.q_4104 +
                a.rather.long.in.fact.very.long.name.q_6101 +
                a.rather.long.in.fact.very.long.name.q_6103 +
                a.rather.long.in.fact.very.long.name.q_7104 +
                a.rather.long.in.fact.very.long.name.q_3109 +
                a.rather.long.in.fact.very.long.name.q_4103 +
                a.rather.long.in.fact.very.long.name.q_2111 +
                a.rather.long.in.fact.very.long.name.q_3107 +
                a.rather.long.in.fact.very.long.name.q_3101 +
                a.rather.long.in.fact.very.long.name.q_3104 +
                a.rather.long.in.fact.very.long.name.q_7107,
                data = a.rather.long.in.fact.very.long.name.for.the.dataframe)
print(summary(a.rather.long.in.fact.very.long.name.for.the.modelA, digits = 2))
plot(a.rather.long.in.fact.very.long.name.for.the.modelA)
plotmo(a.rather.long.in.fact.very.long.name.for.the.modelA)

cat("--a.rather.long.in.fact.very.long.name.for.the...C--\n")
a.rather.long.in.fact.very.long.name.for.the.modelC <-
        earth(x = a.rather.long.in.fact.very.long.name.for.the.dataframe[,-1],
          y = a.rather.long.in.fact.very.long.name.for.the.response,
                  degree = 3)
print(summary(a.rather.long.in.fact.very.long.name.for.the.modelC, digits = 2))
plot(a.rather.long.in.fact.very.long.name.for.the.modelC)
plotmo(a.rather.long.in.fact.very.long.name.for.the.modelC)

a <- earth(survived ~ pclass+sex+age, data=etitanic, degree=2)
print(summary(a))
plotmo(a, caption="plotmo with facs: pclass+sex+age")
plotmo(a, caption="plotmo with facs: pclass+sex+age, all1=T, grid.col=\"gray\"", all1=T, grid.col="gray")
plotmo(a, caption="plotmo with facs: pclass+sex+age, all2=T, col.grid=\"green\"", all2=T, col.grid="green")
plotmo(a, caption="plotmo with facs: pclass+sex+age, all1=T, all2=T, grid=2", all1=T, all2=T, grid.col=2)
plotmo(a, clip=FALSE, degree2=FALSE, caption="plotmo (no degree2) with facs: pclass+sex+age")
plotmo(a, clip=FALSE, grid.levels=list(pclass="2n", sex="ma"),
       caption="plotmo with grid.levels: pclass+sex+age")
# in above tests, all degree2 terms use facs
# now build a model with some degree2 term that use facs, some that don't
a <- earth(survived ~ pclass+age+sibsp, data=etitanic, degree=2)
print(summary(a))
plotmo(a, caption="plotmo with mixed fac and non-fac degree2 terms", persp.border=NA)
plotmo(a, caption="plotmo with mixed fac and non-fac degree2 terms and grid.levels",
       grid.levels=list(pclass="2n", age=20), # test partial matching of grid levels, and numeric preds
       persp.ticktype="d", persp.nticks=2)

# check detection of illegal grid.levels argument
expect.err(try(plotmo(a, grid.levels=list(pcla="1", pclass="2"))), 'illegal grid.levels argument ("pcla" and "pclass" both match "pclass")')
expect.err(try(plotmo(a, grid.levels=list(pclass="1", pcla="2"))), 'illegal grid.levels argument ("pclass" and "pcla" both match "pclass")')
expect.err(try(plotmo(a, grid.levels=list(pcla="nonesuch"))), 'illegal level "nonesuch" for "pclass" in grid.levels (allowed levels are "1st" "2nd" "3rd")')
expect.err(try(plotmo(a, grid.levels=list(pcla="1sx"))), 'illegal level "1sx" for "pclass" in grid.levels (allowed levels are "1st" "2nd" "3rd")')
expect.err(try(plotmo(a, grid.levels=list(pcla=1))), 'illegal level for "pclass" in grid.levels (specify factor levels with a string)')
expect.err(try(plotmo(a, grid.levels=list(pcla=c("ab", "cd")))), "length(pclass) in grid.levels is not 1")
expect.err(try(plotmo(a, grid.levels=list(pcla=NA))), 'pclass in grid.levels is NA')
expect.err(try(plotmo(a, grid.levels=list(pcla=Inf))), 'pclass in grid.levels is infinite')
expect.err(try(plotmo(a, grid.levels=list(pcla=9))), 'illegal level for "pclass" in grid.levels (specify factor levels with a string)')
options(warn=2)
expect.err(try(plotmo(a, grid.levels=list(age="ab"))), 'grid.levels returned class \"character\" for age, so will use the default grid.func for age')
options(warn=1)
expect.err(try(plotmo(a, grid.levels=list(age=NA))), 'age in grid.levels is NA')
expect.err(try(plotmo(a, grid.levels=list(age=Inf))), 'age in grid.levels is infinite')
expect.err(try(plotmo(a, grid.lev=list(age=list(1,2)))), 'length(age) in grid.levels is not 1')

# more-or-less repeat above, but with glm models
a <- earth(survived ~ pclass+age+sibsp, data=etitanic, degree=2, glm=list(family=binomial))
print(summary(a))
plotmo(a, ylim=c(0, 1), caption="plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, ylim=c(0, 1), caption="plotmo glm with mixed fac and non-fac degree2 terms and grid.levels",
       grid.levels=list(pcl="2nd")) # test partial matching of variable name in grid levels
plotmo(a, type="earth", ylim=c(0, 1), caption="type=\"earth\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, type="link", ylim=c(0, 1), clip=FALSE, caption="type=\"link\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, type="class", ylim=c(0, 1), caption="type=\"class\" plotmo glm with mixed fac and non-fac degree2 terms")
plotmo(a, ylim=c(0, 1), caption="default type (\"response\")\nplotmo glm with mixed fac and non-fac degree2 terms")
# now with different type2s
set.seed(2016)
plotmo(a, do.par=FALSE, type2="persp", persp.theta=-20, degree1=FALSE, grid.levels=list(pclass="2nd"))
mtext("different type2s", outer=TRUE, font=2, line=1.5, cex=1)
plotmo(a, do.par=FALSE, type2="contour", degree1=FALSE, grid.levels=list(pclass="2nd"))
plotmo(a, do.par=FALSE, type2="image",   degree1=FALSE, grid.levels=list(pclass="2nd"),
       col.response=as.numeric(etitanic$survived)+2, pt.pch=20)
plotmo(a, do.par=FALSE, type="earth", type2="image", degree1=FALSE,
       grid.levels=list(pclass="2"))

# grid.levels with partdep

set.seed(2018)
x1 <- (1:11) + runif(11)
x2 <- (1:11) + runif(11)
x3 <- as.integer((1:11) + runif(11))
x4 <- runif(11) > .5 # logical
y <- x1 - x2 + x3 + x4
data <- data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4)
lm.x1.x2.x3 <- lm(y ~ x1 + x2 + x3 + x4 + x1*x2 + x1*x3, data=data)
cat("summary(lm.x1.x2.x3):\n")
print(summary(lm.x1.x2.x3))
par(mfrow = c(5, 6), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.5, 0), cex = 0.6, oma=c(0,0,8,0))
plotmo(lm.x1.x2.x3, do.par=0, ylim=c(0,16), pt.col=2,
       caption="row1 default\nrow2 grid.levels=list(x3=15)\nrow3 partdep\nrow4 partdetp grid.levels=list(x3=15)")
plotmo(lm.x1.x2.x3, do.par=0, ylim=c(0,16), pt.col=2, grid.levels=list(x3=15))
plotmo(lm.x1.x2.x3, do.par=0, ylim=c(0,16), pt.col=2, pmethod="partdep")
plotmo(lm.x1.x2.x3, do.par=0, ylim=c(0,16), pt.col=2, pmethod="partdep", grid.levels=list(x3=15))

# check auto type convert in grid.levels
plotmo(lm.x1.x2.x3, degree1="x1", degree2=0, main="x1 (x2=5L))",  ylim=c(0,16), do.par=0, pmethod="partdep", grid.levels=list(x2=15L)) # integer to numeric
plotmo(lm.x1.x2.x3, degree1="x1", degree2=0, main="x1 (x3=5))",   ylim=c(0,16), do.par=0, pmethod="partdep", grid.levels=list(x3=15))  # numeric to integer
plotmo(lm.x1.x2.x3, degree1="x1", degree2=0, main="x1 (x4=1))",   ylim=c(0,16), do.par=0, pmethod="partdep", grid.levels=list(x4=1))   # numeric to logical
expect.err(try(plotmo(lm.x1.x2.x3, degree1="x1", degree2=0, main="x1 (x4=1))",   ylim=c(0,16), do.par=0, pmethod="partdep", grid.levels=list(x4="x"))), "expected a logical value in grid.levels for x4") # char to logical
expect.err(try(plotmo(lm.x1.x2.x3, degree1="x2", do.par=0, pmethod="partdep", grid.levels=list(x1="1"))), "the class \"character\" of \"x1\" in grid.levels does not match its class \"numeric\" in the input data")
par(org.par)

# test vector main

a20 <- earth(O3 ~ humidity + temp + doy, data=ozone1, degree=2, glm=list(family=Gamma))

dopar(2, 2)
plotmo(a20, nrug=-1)

set.seed(2016)
plotmo(a20, nrug=10, caption="Test plotmo with a vector main (and npoints=200)",
       main=c("Humidity", "Temperature", "Day of year", "Humidity: Temperature", "Temperature: Day of Year"),
       col.response="darkgray", pt.pch=".", cex.response=3, npoints=200) # cex.response tests back compat

cat("Expect warning below (missing double titles)\n")
plotmo(a20, nrug=-1, caption="Test plotmo with a vector main (and plain smooth)",
       main=c("Humidity", "Temperature", "Day of year", "Humidity: Temperature", "Temp: Doy"),
       smooth.col="indianred")

cat("Expect warning below (missing single titles)\n")
plotmo(a20, nrug=-1, caption="Test plotmo with a vector main (and smooth args)",
       main=c("Humidity", "Temperature"),
       smooth.col="indianred", smooth.lwd=2, smooth.lty=2, smooth.f=.1,
       col.response="gray", npoints=500)

plotmo(a20, nrug=-1, caption="Test plotmo with pt.pch=paste(1:nrow(ozone1))",
       type2="im",
       col.response=2, pt.cex=.8, pt.pch=paste(1:nrow(ozone1)), npoints=100)

aflip <- earth(O3~vh + wind + humidity + temp, data=ozone1, degree=2)

# test all1 and all2, with and without degree1 and degree2
plotmo(aflip, all2=T, caption="all2=T", npoints=TRUE)
plotmo(aflip, all2=T, degree2=c(4, 2), caption="all2=T, degree2=c(4, 2)")
plotmo(aflip, all1=T, caption="all1=T")
plotmo(aflip, all1=T, degree1=c(3,1), degree2=NA, caption="all1=T, degree1=c(3,1), degree2=NA")

options(warn=2)
expect.err(try(plotmo(aflip, no.such.arg=9)), "(converted from warning) predict.earth ignored argument 'no.such.arg'")
expect.err(try(plotmo(aflip, ycolumn=1)), "(converted from warning) predict.earth ignored argument 'ycolumn'")
expect.err(try(plotmo(aflip, title="abc")), "(converted from warning) predict.earth ignored argument 'title'")
expect.err(try(plotmo(aflip, persp.ticktype="d", persp.ntick=3, tic=3, tick=9)), "(converted from warning) predict.earth ignored argument 'tic'")
expect.err(try(plotmo(aflip, persp.ticktype="d", ntick=3, tic=3)), "(converted from warning) predict.earth ignored argument 'ntick'")
options(warn=1)
# expect.err(try(plotmo(aflip, adj1=8, adj2=9))) # Error : plotmo: illegal argument "adj1"
# expect.err(try(plotmo(aflip, yc=8, x2=9))) # "ycolumn" is no longer legal, use "nresponse" instead
# expect.err(try(plotmo(aflip, persp.ticktype="d", ntick=3, ti=3))) # Error : "title" is illegal, use "caption" instead ("ti" taken to mean "title")
# expect.err(try(plotmo(aflip, persp.ticktype="d", ntick=3, title=3))) # Error : "title" is illegal, use "caption" instead
# expect.err(try(plotmo(aflip, persp.ticktype="d", ntick=3, tit=3, titl=7))) # Error : "title" is illegal, use "caption" instead ("tit" taken to mean "title")
# expect.err(try(plotmo(aflip, zlab="abc"))) # "zlab" is illegal, use "ylab" instead
# expect.err(try(plotmo(aflip, z="abc"))) # "zlab" is illegal, use "ylab" instead ("z" taken to mean "zlab")
expect.err(try(plotmo(aflip, degree1=c(4,1))), "'degree1' is out of range, allowed values are 1 to 2")
# expect.err(try(plotmo(aflip, none.such=TRUE))) # illegal argument "all1"
# expect.err(try(plotmo(aflip, ntick=3, type2="im"))) # the ntick argument is illegal for type2="image"
# expect.err(try(plotmo(aflip, breaks=3, type2="persp"))) # the breaks argument is illegal for type2="persp"
# expect.err(try(plotmo(aflip, breaks=99, type2="cont"))) #  the breaks argument is illegal for type2="contour"

# Test error handling when accessing the original data

lm.bad <- lm.fit(as.matrix(ozone1[,-1]), as.matrix(ozone1[,1]))
expect.err(try(plotmo(lm.bad)), "'lm.bad' is a plain list, not an S3 model")
expect.err(try(plotmo(99)), "'99' is not an S3 model")

x <- matrix(c(1,3,2,4,5,6,7,8,9,10,
              2,3,4,5,6,7,8,9,8,9), ncol=2)

colnames(x) <- c("c1", "c2")
x1 <- x[,1]
x2 <- x[,2]
y <- 3:12
df <- data.frame(y=y, x1=x1, x2=x2)
foo1 <- function()
{
    a.foo1 <- lm(y~x1+x2, model=FALSE)
    x1 <- NULL
    expect.err(try(plotmo(a.foo1)), "cannot get the original model predictors")
}
foo1()
foo2 <- function()
{
    a.foo2 <- lm(y~x1+x2, data=df, model=FALSE)
    df <- 99 # note that df <- NULL here will not cause an error msg
    y <- 99  # also needed else model.frame in plotmo will find the global y
    expect.err(try(plotmo(a.foo2)), "cannot get the original model predictors")
}
foo2()
foo3 <- function()
{
    a.foo3 <- lm(y~x) # lm() builds an lm model for which predict doesn't work
    expect.err(try(plotmo(a.foo3)), "predict returned the wrong length (got 10 but expected 50)")
}
foo3()
foo3a <- function()
{
    a.foo3a <- lm(y~x) # lm() builds an lm model for which predict doesn't work
    # this tests "ngrid1 <- ngrid1 + 1" in plotmo.R
    expect.err(try(plotmo(a.foo3a, ngrid1=nrow(x))), "predict returned the wrong length (got 10 but expected 11)")
}
foo3a()
foo4 <- function()
{
    a.foo4 <- lm(y~x[,1]+x[,2])  # builds an lm model for which predict doesn't work
    # causes 'newdata' had 8 rows but variables found have 10 rows
    expect.err(try(plotmo(a.foo4)), "predict returned the wrong length (got 10 but expected 50)")
}
foo4()
foo5 <- function()
{
    a.foo5 <- lm(y~x1+x2, model=FALSE)
    x1 <- c(1,2,3)
    # causes Error in model.frame.default: variable lengths differ (found for 'x1')
    expect.err(try(plotmo(a.foo5)), "cannot get the original model predictors")
}
foo5()
foo6 <- function()
{
    a.foo6 <- lm(y~x1+x2, model=FALSE)
    y[1] <- NA
    # Error in na.fail.default: missing values in object
    expect.err(try(plotmo(a.foo6, col.response=3)), "cannot get the original model predictors")
}
foo6()
foo7 <- function()
{
    a.foo7 <- lm(y~x1+x2, model=FALSE)
    y[1] <- Inf
    options <- options("warn")
    on.exit(options(warn=options$warn))
    options(warn=2)
    expect.err(try(plotmo(a.foo7, col.response=3)), "non-finite values returned by plotmo_y")
}
foo7()
options(warn=1)
foo8 <- function()
{
    i <- 1
    a.foo8 <- lm(y~x[,i]+x[,2])
    options <- options("warn")
    on.exit(options(warn=options$warn))
    options(warn=2)
    expect.err(try(plotmo(a.foo8)), "Cannot determine which variables to plot in degree2 plots (use all2=TRUE?)")
    options(warn=options$warn)
    expect.err(try(plotmo(a.foo8)), "predict returned the wrong length (got 10 but expected 50)")
}
foo8()
options(warn=1)
foo9 <- function()
{
    my.list <- list(j=2)
    a.foo9 <- lm(y~x[,1]+x[,my.list$j])
    expect.err(try(plotmo(a.foo9)), "cannot get the original model predictors")
}
foo9()
foo9a <- function()
{
    df <- data.frame(y=y, x1=x[,1], x2=x[,2])
    a.foo9a <- lm(y~x1+x2, data=df)
    par(mfrow = c(2, 2), oma=c(0,0,4,0))
    set.seed(2018)
    plotmo(a.foo9a, col.resp=2, do.par=FALSE,
           caption="top two plots should be identical to bottom two plots")
    x2 <- rep(99, length(x2))
    a.foo9b <- lm(y~x1+x2, data=df)
    x2 <- rep(199, length(x2))
    plotmo(a.foo9b, col.resp=2, do.par=FALSE)
}
foo9a()
par(org.par)

foo20.func <- function()
{
    par(mfrow = c(2, 2), oma=c(0,0,4,0))
    foo20 <- lm(y~x1+x2)
    set.seed(2018)
    plotmo(foo20, degree1=1:2, col.resp=2, do.par=FALSE,
           caption="top two plots should be identical to bottom two plots\nbecause we use saved lm$model")
    x1 <- 99
    plotmo(foo20, degree1=1:2, col.resp=2, do.par=FALSE)
}
foo20.func()
par(org.par)

set.seed(1235)
tit <- etitanic
tit <- tit[c(30:80,330:380,630:680), ]
a <- earth(survived~., data=tit, glm=list(family=binomial), degree=2)
plotmo(a, grid.levels=list(sex="ma"),
       caption="smooth: survived, sex=\"m\"    jitter=1",
       smooth.col="indianred", smooth.lwd=2,
       col.response=as.numeric(tit$survived)+2, pt.pch=".", type2="im",
       pt.cex=3, jitter=1) # big jitter
set.seed(1238)
a <- earth(pclass~., data=tit)
plotmo(a, type="class", nresponse=1,
       grid.levels=list(sex="ma"),
       caption="smooth: pclass, sex=\"m\"", SHOWCALL=TRUE,
       smooth.col="indianred", smooth.lwd=2,
       col.response=as.numeric(tit$pclass)+1, type2="im",
       pt.pch=".", pt.cex=3)
plotmo(a, type="class", nresponse=1,
       grid.levels=list(sex="ma"),
       caption="smooth: pclass, sex=\"m\"      jitter=.3", SHOWCALL=TRUE,
       smooth.col="indianred", smooth.lwd=2,
       col.response=as.numeric(tit$pclass)+1, type2="im",
       pt.pch="x", jit=.3) # small jitter
plotmo(a, nresponse=1,
       type="class", grid.levels=list(sex="ma"),
       caption="smooth: pclass, sex=\"m\"",  SHOWCALL=TRUE,
       smooth.col="indianred", smooth.lwd=2,
       col.response=as.numeric(tit$pclass)+1, type2="im",
       pt.pch=paste(1:nrow(tit)))

# test the extend argument

plotmo(a, nresponse=1,             pt.col=2, degree2=0, SHOWCALL=TRUE,
       caption="test extend: extend=0 (reference plot)")
plotmo(a, nresponse=1, extend=.5, pt.col=2, SHOWCALL=TRUE,
       caption="test extend: extend=.5")
plotmo(a, nresponse=1, degree1=0, extend=.2, pt.col=2, SHOWCALL=TRUE) # nothing to plot

a <- earth(survived~pclass+age, data=etitanic, degree=2)
# expect warning: extend=.5 not degree2 plots
plotmo(a, extend=.5, pt.col=2, SHOWCALL=TRUE,
       caption="test extend: extend=.5")

# intercept only models

dopar(2, 2, caption = "intercept-only models")
set.seed(1)
x <- 1:10
y <- runif(length(x))
earth.intercept.only <- earth(x, y)
plotmo(earth.intercept.only, do.par=FALSE, main="earth intercept-only model")
plotmo(earth.intercept.only, do.par=FALSE, col.response=1, pt.pch=20)
# TODO following draws a plot but it shouldn't (very minor bug because int-only model with a bad degree1 spec)
plotmo(earth.intercept.only, do.par=FALSE, degree1=3) # expect warning: 'degree1' specified but no degree1 plots
plotmo(earth.intercept.only, do.par=FALSE, degree1=0) # expect warning: plotmo: nothing to plot
library(rpart)
rpart.intercept.only <- rpart(y~x)
plotmo(rpart.intercept.only, do.par=FALSE, main="rpart.plot intercept-only model")
plotmo(rpart.intercept.only, do.par=FALSE, degree1=0)
par(org.par)

# nrug argument

par(mfrow=c(3,3), mar=c(3,3,3,1), mgp=c(1.5, 0.5, 0))
mod.nrug <- earth(survived~age, data=etitanic)
set.seed(2016)
plotmo(mod.nrug, do.par=0, nrug=-1,   main="nrug=-1")
plotmo(mod.nrug, do.par=0, nrug=TRUE, main="nrug=TRUE")
plotmo(mod.nrug, do.par=0, nrug=10, rug.col=2,            main="nrug=10, rug.col=2")
plotmo(mod.nrug, do.par=0, nrug=5,  rug.col=2, rug.lwd=2, main="nrug=5, rug.col=2, rug.lwd=2")
plotmo(mod.nrug, do.par=0, nrug="density", main="nrug=\"density\"")
plotmo(mod.nrug, do.par=0, nrug="density", density.col=2, density.lwd=2, main="nrug=\"density\"\ndensity.col=2, density.lwd=2")
plotmo(mod.nrug, do.par=0, nrug="density", density.adj=.2, density.col=1, main="nrug=\"density\"\ndensity.adj=.2, density.col=1")
par(org.par)

# a <- earth(ozone1[,3]~ozone1[,1]+ozone1[,2]+ozone1[,4]+ozone1[,5]+ozone1[,6], data=ozone1)
# # TODO fails: actual.nrows=330 expected.nrows=50 fitted.nrows=330
# plotmo(a)

# # TODO following fails in plotmo with
# # Error : get.earth.x from model.matrix.earth from predict.earth: x has 2 columns, expected 4 to match: 1 2 3 Girth
# a <- earth(Volume~poly(Height, degree=3)+Girth, data=trees, subset=4:23, linpreds=TRUE)
# plotmo(a, trace=-1, do.par=FALSE, caption="all three rows should be the same")

source("test.epilog.R")
