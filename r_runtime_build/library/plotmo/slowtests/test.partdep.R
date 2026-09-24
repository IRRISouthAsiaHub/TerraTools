# partdep.test.R: partdep tests for plotmo and plotres

source("test.prolog.R")
library(plotmo)
library(earth)
data(etitanic)

mod <- earth(survived~., data=etitanic, degree=2)

plotmo(mod, caption="plotmo classical")

plotmo(mod, pmethod="partdep", caption="plotmo partdep age")

set.seed(2016)
plotmo(mod, pmethod="apartdep", caption="plotmo apartdep age", do.par=2)

set.seed(2016)
plotmo(mod, pmethod="apartdep", ylim=c(0,1), do.par=0,
       type2="image", pt.col=ifelse(etitanic$survived, "green", "red"),
       degree1=0, degree2=1:3)
par(org.par)

# compare to gbm with an artifical function of variables with a very strong interaction
library(gbm)
n <- 250
set.seed(2016)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- ifelse(x2 > .6, x1-.2, ifelse(x2 > .4, 1 - 1.5 * x1, .3)) + .1 * sin(4 * x3)
data <- data.frame(x1=x1, x2=x2, x3=x3, y=y)
n.trees <- 20
set.seed(2016)
mod <- gbm(y~., data=data, n.trees=n.trees, shrinkage=.1,
           distribution="gaussian", interact=5)
plotmo(mod, degree1=0, persp.ticktype="detailed",
       caption="variables with a strong interaction")
par(mfrow=c(4,4), mar=c(2,3,2,1), mgp=c(1.5, 0.5, 0), oma=c(0,0,6,0))
library(viridis);
image.col <- viridis(100)
ngrid1 <- 50
ngrid2 <- 30
plotmo(mod, pmethod="plot", do.par=0, degree2=2, type2="im", ylim=NULL,
       clip=FALSE, image.col=image.col, ngrid1=ngrid1, ngrid=ngrid2)
title("row1: plotmo classic\nrow2: plotmo apartdep\nrow3: plotmo partdep\nrow4: plot.gbm\n\n\n\n\n\n\n", xpd=NA)
ylim <- c(.21, .40)
set.seed(2016) # for consistent selection of rows for partdep.x
plotmo(mod, pmethod="apartdep",  do.par=0, degree2=2, type2="im", ylim=ylim,
       clip=FALSE, image.col=image.col, ngrid1=ngrid1, ngrid=ngrid2)
plotmo(mod, pmethod="partdep",  do.par=0, degree2=2, type2="im", ylim=ylim,
       clip=FALSE, image.col=image.col, ngrid1=ngrid1, ngrid=ngrid2,
       trace=-1) # check that the pacifier messages are suppressed
plot(mod, i.var=1, n.trees=n.trees, ylim=ylim, continuous.resolution=ngrid1)
plot(mod, i.var=2, n.trees=n.trees, ylim=ylim, continuous.resolution=ngrid1)
plot(mod, i.var=3, n.trees=n.trees, ylim=ylim, continuous.resolution=ngrid1)
# following ignores par(mfrow=c(2,2))
plot(mod, i.var=c(1,3), n.trees=n.trees, continuous.resolution=ngrid2,
     col.regions=image.col, colorkey=FALSE,
     main="gbm plot x1:x3\ncompare to plotmo partdep on previous page")
par(org.par)

#--- compare to gbm and randomForest with a simple regression function

data(scor, package="bootstrap") # some correlated data
n <- 50
x1 <- scale(scor$mec[1:n])
x2 <- scale(scor$vec[1:n])
data <- data.frame(x1=x1, x2=x2)

ngrid1 <- 100

# randomForest, simple regression function
library(randomForest)
data$y <- x1 > -.1 # y depends only on x1 (-.1 hand-tuned to create interesting model surface)
set.seed(2016)
# Expect Warning: The response has five or fewer unique values.  Are you sure you want to do regression?
mod <- randomForest(y~., data=data, ntree=3)
par(mfrow=c(4,2), mar=c(2.5,3,2,1), mgp=c(1.3,0.4,0), oma=c(0,0,7,0))
set.seed(2016) # for consistent jitter of response sites
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       type2="image", main="regression surface",
       pt.col=ifelse(data$y, "green", "red"))
title("RANDOM FOREST SIMPLE REGRESSION MODEL
row1: regression surface
row2: plotmo classic type=response
row3: plotmo partdep type=response
row4: randomForest plot\n\n\n\n\n\n\n",
      xpd=NA, adj=0)
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       persp.border=NA, main="regression surface")
plotmo(mod, pmethod="plotmo",  do.par=0, degree2=0, ngrid1=ngrid1,
       type="response")
plotmo(mod, pmethod="partdep",  do.par=0, degree2=0, ngrid1=ngrid1,
       type="response")
partialPlot(mod, pred.data=data, x.var="x1", n.pt=ngrid1,
            which.class="True")
partialPlot(mod, pred.data=data, x.var="x2", n.pt=ngrid1,
            which.class="True")
par(org.par)

# gbm, simple regression function
library(gbm)
n.trees <- 20
data$y <- x1 > -.6 # y depends only on x1 (-.1 hand-tuned to create interesting model surface)
set.seed(2016)
mod <- gbm(y~., data=data, n.trees=n.trees,
           shrinkage=.1, interaction.depth=4,
           distribution="gaussian")
par(mfrow=c(4,2), mar=c(2.5,3,2,1), mgp=c(1.3,0.4,0), oma=c(0,0,7,0))
set.seed(2016) # for consistent jitter of response sites
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       type2="image", main="regression surface",
       pt.col=ifelse(data$y, "green", "red"))
title("GBM SIMPLE REGRESSION MODEL
row1: regression surface
row2: plotmo classic type=response
row3: plotmo partdep type=response
row4: gbm plot\n\n\n\n\n\n\n",
      xpd=NA, adj=0)
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       persp.border=NA, main="regression surface")
plotmo(mod, pmethod="plotmo",  do.par=0, all1=TRUE, degree2=0,
       ngrid1=ngrid1, type="response")
plotmo(mod, pmethod="partdep",  do.par=0, all1=TRUE, degree2=0,
       ngrid1=ngrid1, type="response")
plot(mod, i.var=1, n.trees=n.trees, continuous.resolution=ngrid1)
plot(mod, i.var=2, n.trees=n.trees, continuous.resolution=ngrid1)
par(org.par)

#--- compare to gbm and randomForest with simple binomial (two class) data

data(scor, package="bootstrap") # some correlated data
n <- 50
x1 <- scale(scor$mec[1:n])
x2 <- scale(scor$vec[1:n])
data <- data.frame(x1=x1, x2=x2)

ngrid1 <- 100

# randomForest, simple binomial (two-class) data
library(randomForest)
# y depends only on x1
# random forest requires a factor for classification (not a logical)
data$y <- factor(as.character(x1 > .4),
                 levels=c("FALSE", "TRUE"),
                 labels=c("False", "True"))
set.seed(2016)
mod <- randomForest(y~., data=data, ntree=3)
par(mfrow=c(4,2), mar=c(2.5,3,2,1), mgp=c(1.3,0.4,0), oma=c(0,0,7,0))
set.seed(2016) # for consistent jitter of response sites
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       type2="image", main="regression surface",
       pt.col=ifelse(data$y=="True", "green", "red"))
title("RANDOM FOREST SIMPLE TWO-CLASS MODEL
row1: regression surface
row2: plotmo partdep type=response (FALSE or TRUE)
row3: plotmo partdep type=prob
row4: randomForest partialPlot (clipped log odds)\n\n\n\n\n\n\n",
      xpd=NA, adj=0)
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       persp.border=NA, main="regression surface")

plotmo(mod, pmethod="partdep",  do.par=0, degree2=0, ngrid1=ngrid1,
       type="response")
plotmo(mod, pmethod="partdep",  do.par=0, degree2=0, ngrid1=ngrid1,
       type="prob", nresponse="True", ylim=c(0,1))
partialPlot(mod, pred.data=data, x.var="x1", n.pt=ngrid1,
            which.class="True", ylim=c(-16,16))
partialPlot(mod, pred.data=data, x.var="x2", n.pt=ngrid1,
            which.class="True", ylim=c(-16,16))
par(org.par)

# gbm, simple binomial (two-class) data
library(gbm)
n.trees <- 10
data$y <- as.numeric(x1 > .6) # y depends only on x1
set.seed(2016)
mod <- gbm(y~., data=data, n.trees=n.trees, shrinkage=.1, interact=4,
           distribution="bernoulli")
par(mfrow=c(4,2), mar=c(2.5,3,2,1), mgp=c(1.3,0.4,0), oma=c(0,0,7,0))
set.seed(2016) # for consistent jitter of response sites
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       type2="image", main="regression surface",
       pt.col=ifelse(data$y, "green", "red"))
title("GBM SIMPLE TWO-CLASS MODEL
row1: regression surface
row2: plotmo partdep type=response (probability)
row4: plotmo partdep type=link (log odds)
row3: gbm plot (log odds)\n\n\n\n\n\n\n",
      xpd=NA, adj=0)
plotmo(mod, degree1=0, ngrid2=100, do.par=0, clip=FALSE,
       persp.border=NA, main="regression surface")
plotmo(mod, pmethod="partdep",  do.par=0, all1=TRUE, degree2=0,
       ngrid1=ngrid1, type="response")
plotmo(mod, pmethod="partdep",  do.par=0, all1=TRUE, degree2=0,
       ngrid1=ngrid1, type="link")
plot(mod, i.var=1, n.trees=n.trees, continuous.resolution=ngrid1)
plot(mod, i.var=2, n.trees=n.trees, continuous.resolution=ngrid1)
par(org.par)

source("test.epilog.R")
