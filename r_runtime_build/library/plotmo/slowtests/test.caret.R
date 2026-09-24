# test.caret.R: test plotmo on caret models
#
# TODO This is a minimal set of tests.

source("test.prolog.R")
library(plotmo)
library(earth)
library(caret)
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
set.seed(2010)
caret.earth.mod <- train(O3~., data=ozone1, method="earth",
                         tuneGrid=data.frame(degree=2, nprune=10))
# SHOWCALL is just a testing thing, so we can see who created the plot on the plot itself
plotmo(caret.earth.mod, trace=1, SHOWCALL=TRUE)
plotmo(caret.earth.mod$finalModel, trace=1, SHOWCALL=TRUE)
plotres(caret.earth.mod, trace=1, SHOWCALL=TRUE)
# plotres(caret.earth.mod$finalModel, trace=1, SHOWCALL=TRUE)

set.seed(2015)
bag <- bagEarth(O3~., data=ozone1, degree=2, B=3)
print(bag$fit)
# pairs are plotted correctly (I think)
plotmo(bag, type="response", trace=1, SHOWCALL=TRUE)
plotres(bag, type="response", trace=1, SHOWCALL=TRUE)

set.seed(2015)
a.bag1 <- bagEarth(trees[,-3], trees[,3], degree=2, B = 3)
plotmo(a.bag1, trace=1, SHOWCALL=TRUE, all2=TRUE, caption="bagEarth, trees")
plotres(a.bag1, trace=1, SHOWCALL=TRUE)

# trace=1 to display "Fixed rank deficient bx by removing 1 term" messages
set.seed(2015)
a.bag3 <- bagEarth(survived~., data=etitanic, degree=2, B=3, trace=1)
plotmo(a.bag3, clip=F, caption="bagEarth, etitanic", trace=1, SHOWCALL=TRUE)
plotres(a.bag3, clip=F, trace=1, SHOWCALL=TRUE)

# following based on example by Max Kuhn on stackoverflow
etit <- etitanic
etit$survived <- factor(ifelse(etit$survived == 1, "yes", "no"),
                       levels = c("yes", "no"))
set.seed(2015)
caret.earth.mod2 <- train(survived ~ .,
            data = etit,
            method = "earth",
            tuneGrid = data.frame(degree = 2, nprune = 9),
            trControl = trainControl(method = "none",
                                     classProbs = TRUE))
# Following gives expected warning (because factors in caret-earth model)
#  Warning: Cannot determine which variables to plot (use all1=TRUE?)
plotmo(caret.earth.mod2, trace=1, SHOWCALL=TRUE)
# changed Sep 2020: following with all2=2 generates the same plot as above (because with warning, above defaults to all2=TRUE)
plotmo(caret.earth.mod2, trace=1, all2=TRUE, SHOWCALL=TRUE, caption="caret.earth.mod2: all2=2")
plotres(caret.earth.mod2, trace=1, SHOWCALL=TRUE)

# Sep 2020: test with a logical variable (check that get.earth.vars.for.plotmo strips "sexTRUE" to "sex")
# following should be exactly the same model as caret.earth.mod2 except for the variable naming for sex
etit.bool <- etitanic
etit.bool$survived <- factor(ifelse(etit.bool$survived == 1, "yes", "no"),
                                levels = c("yes", "no"))
etit.bool$sex <- etit.bool$sex == "male" # to bool
set.seed(2015) # same random seed as above (may not be necessary)
caret.earth.boolfac <- train(survived ~ .,
            data = etit.bool,
            method = "earth",
            tuneGrid = data.frame(degree = 2, nprune = 9),
            trControl = trainControl(method = "none",
                                     classProbs = TRUE))
print(summary(caret.earth.boolfac))

plotmo(caret.earth.boolfac, trace=0, SHOWCALL=TRUE) # Warning: Cannot determine which variables to plot (use all1=TRUE?)
# changed Sep 2020: following with all1=TRUE, all2=TRUE generates the same plot as above
plotmo(caret.earth.boolfac, trace=0, all1=TRUE, all2=TRUE, SHOWCALL=TRUE, caption="caret.earth.mod2: all1=T, all2=T")

data(ozone1)
set.seed(2020)
a <- train(O3 ~ ., data = ozone1,  method = "earth",
            tuneGrid = data.frame(degree = 2, nprune = 14))
plotmo(a, trace=1, SHOWCALL=TRUE)
plotres(a, trace=1, SHOWCALL=TRUE)

cat("=== method=\"svmRadial\" (S4 model wrapped in an S3 model) ===\n")
data(trees)
set.seed(2019)
library(kernlab)
mod <- train(Girth~., data=trees, method="svmRadial",
             trControl=trainControl(method="cv", number=2),
             tuneLength=2, preProcess = c("center", "scale"))
plotres(mod, info=TRUE)
set.seed(2020)
plotmo(mod, pt.col=2, all2=TRUE, pmethod="partdep")

source("test.epilog.R")
