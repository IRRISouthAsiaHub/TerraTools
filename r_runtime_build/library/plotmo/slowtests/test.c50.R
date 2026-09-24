# test.c50.R: c50 tests for plotmo and plotres

source("test.prolog.R")
library(C50)
library(rpart.plot) # for ptitanic, want data with NAs for testing
library(plotmo)
library(earth) # for etitanic
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
tit <- get.tit()

c50.tree.xy <- C5.0(x=tit[,-1], y=tit[,1]) # predict pclass
plotmo(c50.tree.xy, type="prob", nresponse="first", pmethod="apartdep")
plotmo(c50.tree.xy, type="class")
# TODO following gives error: type should be either 'class', 'confidence' or 'prob'
# try(plotmo(c50.tree.xy, type="confidence"))
plotres(c50.tree.xy, type="prob", nresponse="first")

c50.tree.form <- C5.0(pclass~., data=tit) # predict pclass
plotmo(c50.tree.form, type="prob", nresponse="first")
plotmo(c50.tree.form, type="class")
# TODO following gives error: type should be either 'class', 'confidence' or 'prob'
# try(plotmo(c50.tree.form, type="confidence"))
plotres(c50.tree.form, type="prob", nresponse="first")

tit$survived <- factor(ifelse(tit$survived == 1, "yes", "no"),
                       levels = c("yes", "no"))
c50.tree.survived <- C5.0(survived~., data=tit, trials=5) # predict survived
plotmo(c50.tree.survived, type="prob", nresponse="yes")
plotmo(c50.tree.survived, type="class")
# TODO following gives error: type should be either 'class', 'confidence' or 'prob'
# try(plotmo(c50.tree.survived, type="confidence"))
plotres(c50.tree.survived, type="prob", nresponse="yes")

source("test.epilog.R")
