# make.README.figs.R: Build the figures used by plotmo README.md
# The paths below assume that this file is in the plotmo/inst/slowtests directory
# Swindon May 2018

library(plotmo)
library(earth) # for the ozone1 data
data(ozone1)
library(randomForest)
oz <- ozone1[, c("O3", "humidity", "temp")] # small set for illustration
set.seed(2018)
rf.mod <- randomForest(O3 ~ ., data=oz)

# png("../../inst/README-figures/plotmo-randomForest.png", width=460, height=500)
# plotmo(rf.mod, cex.caption=1.5, font.caption=2, oma=c(0,0,5,0),
#        persp.ticktype="detailed", persp.nticks=2)
# dev.off()

# png("../../inst/README-figures/plotres-randomForest.png", width=460, height=530)
# set.seed(2018)
# plotres(rf.mod, cex=1.1, cex.caption=1.5, font.caption=2, oma=c(1,0,3,0))
# dev.off()

# png("../../inst/README-figures/plotres-glmnet-gbm.png", width=700, height=400)
# par(mfrow=c(1,2), oma=c(1,0,0,0))
# library(glmnet)
# set.seed(2016)
# x <- matrix(rnorm(100 * 10), 100, 10)
# y <- x[,1] + x[,2] + 3 * rnorm(100) # y depends only on x[,1] and x[,2]
# mod <- glmnet(x, y)
# plotres(mod, which=1, predict.s=0.25, cex=1.2, pt.cex=.8)
# title("glmnet model\n\n\n")

# library(gbm)
# library(earth); data(ozone1) # get the ozone data
# set.seed(2017)
# oz <- ozone1[sample.int(n=nrow(ozone1)),] # randomize row order for train.fraction
# gbm.mod <- gbm(O3~., data=oz, distribution="gaussian", interaction.depth=2,
#                shrinkage=.01, train.fraction=.8, cv.folds=10, n.trees=3000)
# plotres(gbm.mod, which=1)
# title("gbm model\n\n", xpd=NA)
# dev.off()
