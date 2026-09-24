# test.linmod.R: test example S3 model at http://www.milbo.org/doc/linmod.R

source("test.prolog.R")
source("linmod.R")         # linear model code (http://www.milbo.org/doc/linmod.R)
source("linmod.methods.R") # additional method functions for linmod
options(warn=1) # print warnings as they occur

almost.equal <- function(x, y, max=1e-8)
{
    stopifnot(max >= 0 && max < .01)
    length(x) == length(y) && max(abs(x - y)) < max
}
# check that linmod model matches reference lm model in all essential details
check.lm <- function(fit, ref, newdata=trees[3:5,],
                     check.coef.names=TRUE,
                     check.casenames=TRUE,
                     check.newdata=TRUE,
                     check.sigma=TRUE)
{
    check.names <- function(fit.names, ref.names)
    {
        if(check.casenames &&
        # lm always adds rownames even if "1", "2", "3": this seems
        # wasteful and not particulary helpful, so linmod doesn't do
        # this, hence the first !isTRUE(all.equal) below
           !isTRUE(all.equal(ref.names, paste(1:length(ref.names)))) &&
           !isTRUE(all.equal(fit.names, ref.names))) {
            print(fit.names)
            print(ref.names)
            stop(deparse(substitute(fit.names)), " != ",
                 deparse(substitute(ref.names)))
        }
    }
    cat0("check ", deparse(substitute(fit)), " vs ",
         deparse(substitute(ref)), "\n")

    stopifnot(coef(fit) == coef(ref))
    if(check.coef.names)
        stopifnot(identical(names(coef(fit)), names(coef(ref))))

    stopifnot(identical(dim(fit$coefficients), dim(ref$coefficients)))
    stopifnot(length(fit$coefficients) == length(ref$coefficients))
    stopifnot(almost.equal(fit$coefficients, ref$coefficients))

    stopifnot(identical(dim(fit$residuals), dim(ref$residuals)))
    stopifnot(length(fit$residuals) == length(ref$residuals))
    stopifnot(almost.equal(fit$residuals, ref$residuals))

    stopifnot(identical(dim(fit$fitted.values), dim(ref$fitted.values)))
    stopifnot(length(fit$fitted.values) == length(ref$fitted.values))
    stopifnot(almost.equal(fit$fitted.values, ref$fitted.values))

    stopifnot(identical(fit$rank, ref$rank))

    if(!is.null(fit$vcov) && !is.null(ref$vcov)) {
        stopifnot(identical(dim(fit$vcov), dim(ref$vcov)))
        stopifnot(length(fit$vcov) == length(ref$vcov))
        stopifnot(almost.equal(fit$vcov, ref$vcov))
    }
    if(check.sigma) {
        ref.sigma <- ref$sigma
        if(is.null(ref.sigma)) # in lm models, sigma is only available from summary()
            ref.sigma <- summary(ref)$sigma
        stopifnot(almost.equal(fit$sigma, ref.sigma))
    }
    stopifnot(almost.equal(fit$df.residual, ref$df.residual))

    stopifnot(almost.equal(fitted(fit), fitted(ref)))
    check.names(names(fitted(fit)), names(fitted(ref)))

    stopifnot(almost.equal(residuals(fit), residuals(ref)))
    check.names(names(residuals(fit)), names(residuals(ref)))

    stopifnot(almost.equal(predict(fit), predict(ref)))
    check.names(names(predict(fit)), names(predict(ref)))
    if(check.newdata) {
        stopifnot(almost.equal(predict(fit, newdata=newdata),
                               predict(ref, newdata=newdata)))
        check.names(names(predict(fit, newdata=newdata)),
                    names(predict(ref, newdata=newdata)))
    }
}
tr <- trees # trees data but with rownames
rownames(tr) <- paste("tree", 1:nrow(trees), sep="")

linmod.form.Volume.tr <- linmod(Volume~., data=tr)
cat0("==print(summary(linmod.form.Volume.tr))\n")
print(summary(linmod.form.Volume.tr))
lm.Volume.tr <- lm(Volume~., data=tr)
check.lm(linmod.form.Volume.tr, lm.Volume.tr)
stopifnot(almost.equal(predict(linmod.form.Volume.tr, newdata=data.frame(Girth=10, Height=80)),
                       16.234045, max=1e-5))
stopifnot(almost.equal(predict(linmod.form.Volume.tr, newdata=as.matrix(tr[1:3,])),
                       c(4.8376597, 4.5538516, 4.8169813), max=1e-5))
# character new data (instead of numeric)
newdata.allchar <- as.data.frame(matrix("blank", ncol=3, nrow=3))
colnames(newdata.allchar) <- colnames(trees)
expect.err(try(predict(lm.Volume.tr, newdata=newdata.allchar)),
               "variables 'Girth', 'Height' were specified with different types from the fit")
expect.err(try(predict(linmod.form.Volume.tr, newdata=newdata.allchar)),
              "variables 'Girth', 'Height' were specified with different types from the fit")

linmod.xy.Volume.tr <- linmod(tr[,1:2], tr[,3,drop=FALSE])                         # x=data.frame y=data.frame
cat0("==print(summary(linmod.xy.Volume.tr))\n")
print(summary(linmod.xy.Volume.tr))
newdata.2col <- trees[3:5,1:2]
check.lm(linmod.xy.Volume.tr, lm.Volume.tr, newdata=newdata.2col)
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=data.frame(Girth=10, Height=80)),
                       16.234045, max=1e-5))
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=tr[1:3,1:2]),
                       c(4.8376597, 4.5538516, 4.8169813), max=1e-5))

linmod50.xy.Volume.tr <- linmod(as.matrix(tr[,1:2]), as.matrix(tr[,3,drop=FALSE])) # x=matrix y=matrix
check.lm(linmod50.xy.Volume.tr, lm.Volume.tr, newdata=newdata.2col)
linmod51.xy.Volume.tr <- linmod(tr[,1:2], tr[,3])                                  # x=data.frame y=vector
check.lm(linmod51.xy.Volume.tr, lm.Volume.tr, newdata=newdata.2col)
linmod52.xy.Volume.tr <- linmod(as.matrix(tr[,1:2]), tr[,3])                       # x=matrix y=vector
check.lm(linmod52.xy.Volume.tr, lm.Volume.tr, newdata=newdata.2col)

# newdata can be a vector
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=c(8.3, 70)),
                       4.8376597, max=1e-5))
stopifnot(almost.equal(predict(linmod.xy.Volume.tr,
                       newdata=c(8.3, 8.6, 70, 65)), # 4 element vector, byrow=FALSE
                       c(4.8376597, 4.5538516), max=1e-5))
options(warn=1) # print warnings as they occur
# expect Warning: data length [3] is not a sub-multiple or multiple of the number of rows [2]
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=c(8.3, 9, 70)), # 3 element vector
                      c(4.8376597, -12.7984291), max=1e-5))
options(warn=2) # treat warnings as errors

stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=as.matrix(data.frame(Girth=10, Height=80))),
                       16.234045, max=1e-5))
# column names in newdata are ignored for linmod.default models
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=data.frame(name1.not.in.orig.data=10, name2.not.in.orig.datax2=80)),
                       16.234045, max=1e-5))
# note name reversed below but names still ignored, same predict result as c(Girth=10, Height=80)
stopifnot(almost.equal(predict(linmod.xy.Volume.tr, newdata=data.frame(Height=10, Girth=80)),
                       16.234045, max=1e-5))

cat0("==print.default(linmod.form.Volume.tr)\n")
print.default(linmod.form.Volume.tr)

cat0("==check single x variable\n")
linmod1a.form <- linmod(Volume~Height, data=tr)
cat0("==print(summary(linmod1a.form))\n")
print(summary(linmod1a.form))
lma.tr <- lm(Volume~Height, data=tr)
check.lm(linmod1a.form, lma.tr)

stopifnot(almost.equal(predict(linmod1a.form, newdata=data.frame(Height=80)),
                       36.34437, max=1e-5))
stopifnot(almost.equal(predict(linmod1a.form, newdata=data.frame(Girth=99, Height=80)),
                       36.34437, max=1e-5))
stopifnot(almost.equal(predict(linmod1a.form, newdata=as.matrix(tr[1:3,])),
                       c(20.91087, 13.19412, 10.10742), max=1e-5))

linmod1a.xy <- linmod(tr[,2,drop=FALSE], tr[,3])
cat0("==print(summary(linmod1a.xy))\n")
print(summary(linmod1a.xy))
check.lm(linmod1a.xy, lma.tr, newdata=trees[3:5,2,drop=FALSE])
check.lm(linmod1a.xy, lma.tr, newdata=trees[3:5,2,drop=TRUE],
         check.newdata=FALSE) # needed because predict.lm gives 'data' must be a data.frame, environment, or list
stopifnot(almost.equal(predict(linmod1a.xy, newdata=trees[3:5,2,drop=FALSE]),
                       predict(linmod1a.xy, newdata=trees[3:5,2,drop=TRUE])))
stopifnot(almost.equal(predict(linmod1a.xy, newdata=data.frame(Height=80)),
                       36.34437, max=1e-5))
stopifnot(almost.equal(predict(linmod1a.xy, newdata=tr[1:3,2]),
                       c(20.91087, 13.19412, 10.10742), max=1e-5))
stopifnot(almost.equal(predict(linmod1a.xy, newdata=as.matrix(data.frame(Height=80))),
                       36.34437, max=1e-5))

# check that extra fields in predict newdata are ok with formula models
stopifnot(almost.equal(predict(linmod.form.Volume.tr, newdata=data.frame(Girth=10, Height=80, extra=99)),
                       predict(lm.Volume.tr,          newdata=data.frame(Girth=10, Height=80))))
stopifnot(almost.equal(predict(linmod.form.Volume.tr, newdata=data.frame(Girth=10, Height=80, extra=99)),
                       predict(lm.Volume.tr,          newdata=data.frame(Girth=10, Height=80, extra=99))))
# check that extra fields in predict newdata are not ok with x,y models
expect.err(try(predict(linmod.xy.Volume.tr, newdata=data.frame(Girth=10, Height=80, extra=99))),
              "ncol(newdata) is 3 but should be 2")

# missing variables in newdata
expect.err(try(predict(linmod.form.Volume.tr, newdata=data.frame(Girth=10))),
               "object 'Height' not found")
expect.err(try(predict(linmod.form.Volume.tr, newdata=c(8.3, 70))),
               "object 'Girth' not found")
expect.err(try(predict(lm.Volume.tr,          newdata=data.frame(Girth=10))),
               "object 'Height' not found")
expect.err(try(predict(linmod.xy.Volume.tr,            newdata=data.frame(Girth=10))),
               "ncol(newdata) is 1 but should be 2")

# check that rownames got propagated
stopifnot(names(linmod.form.Volume.tr$residuals)[1] == "tree1")
stopifnot(names(linmod.form.Volume.tr$fitted.values)[3] == "tree3")
stopifnot(names(linmod.xy.Volume.tr$residuals)[1] == "tree1")
stopifnot(names(linmod.xy.Volume.tr$fitted.values)[3] == "tree3")
stopifnot(!is.null(names(linmod.xy.Volume.tr$residuals)))
stopifnot(!is.null(names(linmod.xy.Volume.tr$fitted.values)))
cat0("==print.default(linmod.xy.Volume.tr)\n")
print.default(linmod.xy.Volume.tr)

# check that we don't artificially add rownames when no original rownames
linmod1a.xy <- linmod(trees[,1:2], trees[,3])
stopifnot(is.null(names(linmod1a.xy$residuals)))
stopifnot(is.null(names(linmod1a.xy$fitted.values)))

cat0("==example plots\n")

library(plotmo)
data(trees)

linmod.form.Volume.trees <- linmod(Volume~., data=trees)
print(linmod.form.Volume.trees)
print(summary(linmod.form.Volume.trees))

linmod1.xy <- linmod(trees[,1:2], trees[,3])
print(linmod1.xy)
print(summary(linmod1.xy))

plotmo(linmod.form.Volume.trees)
plotmo(linmod1.xy)

plotres(linmod.form.Volume.trees)
plotres(linmod1.xy)

cat0("==test keep arg\n")

trees1 <- trees
linmod.form.Volume.trees.keep <- linmod(Volume~., data=trees1, keep=TRUE)
print(summary(linmod.form.Volume.trees.keep))
print(head(linmod.form.Volume.trees.keep$data))
stopifnot(dim(linmod.form.Volume.trees.keep$data) == c(nrow(trees1), ncol(trees1)))
trees1 <- NULL # destroy orginal data so plotmo has to use keep data
plotmo(linmod.form.Volume.trees.keep, pt.col=3)
plotres(linmod.form.Volume.trees.keep)

linmod.xy.keep <- linmod(trees[,1:2], trees[,3], keep=TRUE)
print(summary(linmod.xy.keep))
print(head(linmod.xy.keep$x))
stopifnot(dim(linmod.xy.keep$x) == c(nrow(trees), 2))
stopifnot(class(linmod.xy.keep$x)[1] == "matrix")
print(head(linmod.xy.keep$y))
stopifnot(dim(linmod.xy.keep$y) == c(nrow(trees), 1))
stopifnot(class(linmod.xy.keep$y)[1] == "matrix")
linmod.xy.keep$call <- NULL # trick to force use of x and y in plotmo
plotmo(linmod.xy.keep, pt.col=3)
plotres(linmod.xy.keep)

check.lm(linmod.form.Volume.trees.keep, linmod.xy.keep, check.casenames=FALSE, check.newdata=FALSE)

cat0("==test keep arg with vector x\n")

n <- 20
linmod.vecx.form.keep <- linmod(Volume~Height, data=trees[1:n,], keep=TRUE)
print(summary(linmod.vecx.form.keep))
print(head(linmod.vecx.form.keep$data))
stopifnot(dim(linmod.vecx.form.keep$data) == c(n, ncol(trees)))
stopifnot(class(linmod.vecx.form.keep$data) == class(trees))
plotmo(linmod.vecx.form.keep, pt.col=3)
plotres(linmod.vecx.form.keep)

linmod.vecx.xy.keep <- linmod(trees[1:n,2], trees[1:n,3], keep=TRUE)
print(summary(linmod.vecx.xy.keep))
print(head(linmod.vecx.xy.keep$x))
stopifnot(dim(linmod.vecx.xy.keep$x) == c(n, 1))
stopifnot(class(linmod.vecx.xy.keep$x)[1] == "matrix")
print(head(linmod.vecx.xy.keep$y))
stopifnot(dim(linmod.vecx.xy.keep$y) == c(n, 1))
stopifnot(class(linmod.vecx.xy.keep$y)[1] == "matrix")
linmod.vecx.xy.keep$call <- NULL # trick to force use of x and y in plotmo
plotmo(linmod.vecx.xy.keep, pt.col=3)
plotres(linmod.vecx.xy.keep)

check.lm(linmod.vecx.form.keep, linmod.vecx.xy.keep, newdata=trees[3:5,2,drop=FALSE],
         check.coef.names=FALSE, check.casenames=FALSE)

cat0("==test model building with assorted numeric args\n")

x <- tr[,1:2]
y <- tr[,3]
cat0("class(x)=", class(x), " class(y)=", class(y), "\n") # class(x)=data.frame class(y)=numeric
linmod2.xy <- linmod(x, y)
check.lm(linmod2.xy, lm.Volume.tr, newdata=newdata.2col)

# check consistency with lm
expect.err(try(linmod(y~x)), "invalid type (list) for variable 'x'")
expect.err(try(lm(y~x)),     "invalid type (list) for variable 'x'")

linmod3.xy <- linmod(as.matrix(x), as.matrix(y))
check.lm(linmod3.xy, lm.Volume.tr, newdata=newdata.2col)

linmod4.form <- linmod(y ~ as.matrix(x))
lm4 <- lm(y ~ as.matrix(x))
check.lm(linmod4.form, lm4, check.newdata=FALSE)
stopifnot(coef(linmod4.form)  == coef(lm.Volume.tr),
          gsub("as.matrix(x)", "", names(coef(linmod4.form)), fixed=TRUE)  == names(coef(lm.Volume.tr)))

xm <- as.matrix(x)
cat0("class(xm)=", class(xm), " class(y)=", class(y), "\n") # class(xm)=matrix class(y)=numeric
linmod5.form <- linmod(y ~ xm)
lm5 <- lm(y ~ xm)
check.lm(linmod5.form, lm5, check.newdata=FALSE)
stopifnot(coef(linmod5.form)  == coef(lm.Volume.tr),
          gsub("xm", "", names(coef(linmod5.form)), fixed=TRUE)  == names(coef(lm.Volume.tr)))

cat0("==test correct use of global x1 and y1, and of predict error handling\n")
x1 <- tr[,1]
y1 <- tr[,3]
cat0("class(x1)=", class(x1), " class(y1)=", class(y1), "\n") # class(x1)=numeric class(y1)=numeric
linmod.y1.x1 <- linmod(y1~x1)
lm1 <- lm(y1~x1)
linmod6.xy <- linmod(x1, y1)

newdata.x1 <- trees[3:5,1,drop=FALSE]
colnames(newdata.x1) <- "x1"
stopifnot(almost.equal(predict(linmod.y1.x1, newdata=newdata.x1),
          c(7.63607739644657, 16.24803331528098, 17.26120459984973)))

check.lm(linmod6.xy, linmod.y1.x1, newdata=x1[3:5],
         check.newdata=FALSE, # TODO needed because linmod.y1.x1 ignores newdata(!)
         check.coef.names=FALSE, check.casenames=FALSE)
print(predict(linmod6.xy, newdata=x1[3:5]))
stopifnot(almost.equal(predict(linmod6.xy, newdata=x1[3]), 7.63607739644657))

stopifnot(coef(linmod6.xy) == coef(linmod.y1.x1)) # names(coef(linmod.y1.x1) are "(Intercept)" "x1"
stopifnot(names(coef(linmod6.xy)) == c("(Intercept)", "V1"))

# following checks some confusing behaviour of predict.lm
options(warn=2) # treat warnings as errors
expect.err(try(predict(lm1,    newdata=trees[3:5,1,drop=FALSE])),
           "'newdata' had 3 rows but variables found have 31 rows")
expect.err(try(predict(lm1,    newdata=trees[3:5,1,drop=TRUE])),
           "'data' must be a data.frame, environment, or list")

# following checks messages when missing variables in newdata

options(warn=2) # treat warnings as errors to check that we get a warning in stats::model.frame
expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=FALSE])),
           "(converted from warning) 'newdata' had 3 rows but variables found have 31 rows")
expect.err(try(predict(lm1, newdata=trees[3:5,1,drop=FALSE])),
           "(converted from warning) 'newdata' had 3 rows but variables found have 31 rows")
expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=TRUE])),
           "(converted from warning) 'newdata' had 3 rows but variables found have 31 rows")

# following checks predict.linmod error messages when missing variables
# (it tries to give better error messages than predict.lm)

options(warn=1) # print warnings as they occur to test stop() in linmod.R::process.newdata.formula
expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=FALSE])),
           "newdata has 3 rows but model.frame returned 31 rows (variable 'x1' may be missing from newdata)")
expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=TRUE])),
           "newdata has 3 rows but model.frame returned 31 rows (variable 'x1' may be missing from newdata)")
options(warn=2) # back to treating warnings as errors

# test old version of linmod.R (pre Sep 2020)
#
# expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=FALSE])),
#            "variable 'x1' is missing from newdata")
# expect.err(try(predict(lm1, newdata=trees[3:5,1,drop=FALSE])),
#            "(converted from warning) 'newdata' had 3 rows but variables found have 31 rows")
# expect.err(try(predict(linmod.y1.x1, newdata=trees[3:5,1,drop=TRUE])),
#            "variable 'x1' is missing from newdata")

linmod6.form <- linmod(y1~x1)
check.lm(linmod6.form, linmod.y1.x1, check.newdata=FALSE)

newdata <- trees[5:6,]
colnames(newdata) <- c("Girth", "Height", "Volume999") # doesn't matter what we call the response
stopifnot(identical(predict(linmod.form.Volume.tr, newdata=newdata),
                    predict(linmod.form.Volume.tr, newdata=trees[5:6,])))
newdata <- trees[5:6,3:1] # reverse columns and their colnames
colnames(newdata) <- c("Volume", "Height", "Girth")
stopifnot(identical(predict(linmod.form.Volume.tr, newdata=newdata),
                    predict(linmod.form.Volume.tr, newdata=trees[5:6,])))
newdata <- trees[5:6,2:1] # reverse columns and their colnames, delete response column
colnames(newdata) <- c("Height", "Girth")
stopifnot(identical(predict(linmod.form.Volume.tr, newdata=newdata),
                    predict(linmod.form.Volume.tr, newdata=trees[5:6,])))
stopifnot(identical(predict(linmod.form.Volume.tr, newdata=as.matrix(trees[5:6,])), # allow matrix newdata
                    predict(linmod.form.Volume.tr, newdata=trees[5:6,])))
newdata <- trees[5:6,]
colnames(newdata) <- c("Girth99", "Height", "Volume")
expect.err(try(predict(linmod.form.Volume.tr, newdata=newdata)),
           "object 'Girth' not found")
colnames(newdata) <- c("Girth", "Height99", "Volume")
expect.err(try(predict(linmod.form.Volume.tr, newdata=newdata)),
           "object 'Height' not found")

cat0("==check integer input (sibsp is an integer)\n")

library(earth) # for etitanic data
data(etitanic)
tit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
tit$survived <- tit$survived != 0 # convert to logical
rownames(tit) <- paste("pas", 1:nrow(tit), sep="")
cat0(paste(colnames(tit), "=", sapply(tit, class), sep="", collapse=", "), "\n")

linmod7.xy <- linmod(tit$age, tit$sibsp)
lm7 <- lm.fit(cbind(1, tit$age), tit$sibsp)
stopifnot(coef(linmod7.xy) == coef(lm7)) # coef names will differ

linmod7.form <- linmod(sibsp~age, data=tit)
lm7.form  <- lm(sibsp~age, data=tit)
check.lm(linmod7.form, lm7.form, newdata=tit[3:5,])

linmod8.xy <- linmod(tit$sibsp, tit$age)
lm8 <- lm.fit(cbind(1, tit$sibsp), tit$age)
stopifnot(coef(linmod8.xy) == coef(lm8)) # coef names will differ

linmod8.form <- linmod(age~sibsp, data=tit)
lm8.form  <- lm(age~sibsp, data=tit)
check.lm(linmod8.form, lm8.form, newdata=tit[3:5,])

# drop=FALSE so response is a data frame
linmod1a.xy <- linmod(trees[,1:2], trees[, 3, drop=FALSE])
print(linmod1a.xy)
print(summary(linmod1a.xy))
plotres(linmod1a.xy) # plot caption shows response name "Volume"

cat0("==test model building with assorted non-numeric args\n")

library(earth) # for etitanic data
data(etitanic)
etit <- etitanic[seq(1, nrow(etitanic), by=60), ] # small set of data for tests (18 cases)
etit$survived <- etit$survived != 0 # convert to logical
rownames(etit) <- paste("pas", 1:nrow(etit), sep="")
cat0(paste(colnames(etit), "=", sapply(etit, class), sep="", collapse=", "), "\n")

lm9 <- lm(survived~., data=etit)
linmod9.form <- linmod(survived~., data=etit)
check.lm(linmod9.form, lm9, newdata=etit[3:5,])

# change class of pclass to numeric
etit.pclass.numeric <- etit
etit.pclass.numeric$pclass <- as.numeric(etit$pclass)
expect.err(try(predict(lm9,       newdata=etit.pclass.numeric)),
              "(converted from warning) variable 'pclass' is not a factor")
expect.err(try(predict(linmod9.form, newdata=etit.pclass.numeric)),
              "(converted from warning) variable 'pclass' is not a factor")

# change class of age to factor
etit.age.factor <- etit
etit.age.factor$age <- etit$pclass
expect.err(try(predict(lm9,       newdata=etit.age.factor)),
              "variable 'age' was fitted with type \"numeric\" but type \"factor\" was supplied")
expect.err(try(predict(linmod9.form, newdata=etit.age.factor)),
              "variable 'age' was fitted with type \"numeric\" but type \"factor\" was supplied")

# predict for formula model ignores extra column(s) in newdata
etit.extra.col <- etit
etit.extra.col$extra <- etit$sibsp
stopifnot(identical(predict(lm9, newdata=etit), predict(lm9, newdata=etit.extra.col)))
stopifnot(identical(predict(linmod9.form, newdata=etit), predict(linmod9.form, newdata=etit.extra.col)))
etit.extra.col$extra2 <- etit$sibsp
stopifnot(identical(predict(lm9, newdata=etit), predict(lm9, newdata=etit.extra.col)))
stopifnot(identical(predict(linmod9.form, newdata=etit), predict(linmod9.form, newdata=etit.extra.col)))

# predict for formula model doesn't care if columns in different order
etit.different.col.order <- etit[,ncol(etit):1] # reverse order of columns
stopifnot(identical(predict(lm9, newdata=etit), predict(lm9, newdata=etit.different.col.order)))
stopifnot(identical(predict(linmod9.form, newdata=etit), predict(linmod9.form, newdata=etit.different.col.order)))

# linmod.default, non numeric x (factors in x)
expect.err(try(linmod(etit[c(1,3,4,5,6)], etit[,"survived"])),
              "non-numeric column in 'x'")
expect.err(try(linmod.fit(etit[c(1,3,4,5,6)], etit[,"survived"])),
              "'x' is not a matrix or could not be converted to a matrix")
# lousy error message from lm.fit
expect.err(try(lm.fit(etit[,c(1,3,4,5,6)], etit[,"survived"])),
              "INTEGER() can only be applied to a 'integer', not a 'NULL'")

expect.err(try(linmod(data.matrix(cbind("(Intercept)"=1, etit[,c(1,3,4,5,6)])), etit[,"survived"])),
              "column name \"(Intercept)\" in 'x' is duplicated")
linmod9a.xy <- linmod(data.matrix(etit[,c(1,3,4,5,6)]), etit[,"survived"])
lm9.fit <- lm.fit(data.matrix(cbind("(Intercept)"=1, etit[,c(1,3,4,5,6)])), etit[,"survived"])
stopifnot(coef(linmod9a.xy) == coef(lm9.fit))
stopifnot(names(coef(linmod9a.xy)) == names(coef(lm9.fit)))
expect.err(try(predict(linmod9a.xy, newdata=etit.age.factor[,c(1,3,4,5,6)])), "non-numeric column in 'newdata'")
expect.err(try(predict(linmod9a.xy, newdata=etit[,c(1,3,4,5)])), "ncol(newdata) is 4 but should be 5")
expect.err(try(predict(linmod9a.xy, newdata=etit[,c(1,3,4,5,6,6)])), "ncol(newdata) is 6 but should be 5")

# linmod.formula, logical response
data.logical.response <- data.frame(etit[1:6,c("age","sibsp","parch")], response=c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
linmod9b.form <- linmod(response~., data=data.logical.response)
print(linmod9b.form)
lm9b.form <- lm(response~., data=data.logical.response)
check.lm(linmod9b.form, lm9b.form, newdata=data.logical.response[2,,drop=FALSE])

# linmod.formula, factor response (not allowed)
data.fac.response <- data.frame(etit[1:6,c("age","sibsp","parch")], response=factor(c("a", "a", "b", "a", "b", "b")))
expect.err(try(linmod(response~., data=data.fac.response)), "'y' is not numeric or logical")
# lm.formula
expect.err(try(lm(response~., data=data.fac.response)),
              "(converted from warning) using type = \"numeric\" with a factor response will be ignored")

# linmod.formula, string response (not allowed)
data.string.response <- data.frame(etit[1:6,c("age","sibsp","parch")], response=c("a", "a", "b", "a", "b", "b"))
expect.err(try(linmod(response~., data=data.string.response)), "'y' is not numeric or logical")
# lm.formula
expect.err(try(lm(response~., data=data.string.response)),
              "(converted from warning) NAs introduced by coercion")

# linmod.default, logical response
linmod9b.xy <- linmod(etit[1:6,c("age","sibsp","parch")], c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
print(linmod9b.xy)
# lm.fit, logical response (lousy error message from lm.fit)
expect.err(try(lm.fit(etit[1:6,c("age","sibsp","parch")], c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))),
              "INTEGER() can only be applied to a 'integer', not a 'NULL'")
# linmod.default, factor response
expect.err(try(linmod(etit[1:6,c("age","sibsp","parch")], factor(c("a",
              "a", "b", "a", "b", "b")))), "'y' is not numeric or logical")
# linmod.default, string response
expect.err(try(linmod(etit[1:6,c("age","sibsp","parch")], c("a",
              "a", "b", "a", "b", "b"))), "'y' is not numeric or logical")
# lm.fit, string and factor responses (lousy error messages from lm.fit)
expect.err(try(lm.fit(etit[1:6,c("age","sibsp","parch")], factor(c("a",
              "a", "b", "a", "b", "b")))), "INTEGER() can only be applied to a 'integer', not a 'NULL'")
expect.err(try(lm.fit(etit[1:6,c("age","sibsp","parch")], c("a",
              "a", "b", "a", "b", "b"))), "INTEGER() can only be applied to a 'integer', not a 'NULL'")

options(warn=2) # treat warnings as errors
expect.err(try(lm(pclass~., data=etit)), "using type = \"numeric\" with a factor response will be ignored")
expect.err(try(linmod(pclass~., data=etit)), "'y' is not numeric or logical")

options(warn=1) # print warnings as they occur
lm10 <- lm(pclass~., data=etit) # will give warnings
options(warn=2) # treat warnings as errors
linmod10.form <- linmod(as.numeric(pclass)~., data=etit)
stopifnot(coef(linmod10.form) == coef(lm10))
stopifnot(names(coef(linmod10.form)) == names(coef(lm10)))
# check.lm(linmod10.form, lm10) # fails because lm10 fitted is all NA

expect.err(try(linmod(pclass~., data=etit)), "'y' is not numeric or logical")
expect.err(try(linmod(etit[,-1], etit[,1])), "non-numeric column in 'x'")
expect.err(try(linmod(1:10, paste(1:10))), "'y' is not numeric or logical")

linmod10a.form <- linmod(survived~pclass, data=etit)
lm10a <- lm(survived~pclass, data=etit)
check.lm(linmod10a.form, lm10a, newdata=etit[3:5,])

expect.err(try(linmod(etit[,"pclass"], etit[,"age"])), "non-numeric column in 'x'")

expect.err(try(linmod(paste(1:10), 1:10)), "non-numeric column in 'x'")

lm11 <- lm(as.numeric(pclass)~., data=etit)
linmod11.form <- linmod(as.numeric(pclass)~., data=etit)
check.lm(linmod11.form, lm11, newdata=etit[3:5,])

# logical data (not numeric)
bool.data <- data.frame(x=rep(c(TRUE, FALSE, TRUE), length.out=10),
                        y=rep(c(TRUE, FALSE, FALSE), length.out=10))
lm12 <- lm(y~x, data=bool.data)
linmod12.form <- linmod(y~x, data=bool.data)
check.lm(linmod12.form, lm12, newdata=bool.data[3:5,1],
         check.newdata=FALSE) # needed because predict.lm gives invalid type (list) for variable 'x'
linmod12.xy <- linmod(bool.data$x, bool.data$y)
# hack: delete mismatching names so check.lm() doesn't fail
names(lm12$coefficients) <- NULL     # were "(Intercept)" "xTRUE"
names(linmod12.xy$coefficients) <- NULL # were "(Intercept)" "V1"
check.lm(linmod12.xy, lm12, newdata=bool.data[3:5,1],
         check.newdata=FALSE, # needed because predict.lm gives invalid 'envir' argument of type 'logical'
         check.casenames=FALSE)

cat0("==check use of functions in arguments to linmod\n")

identfunc <- function(x) x
lm10 <- lm(        sqrt(survived) ~ I(age^2) + as.numeric(sex), data=identfunc(etit))
linmod10 <- linmod(sqrt(survived) ~ I(age^2) + as.numeric(sex), data=identfunc(etit))
print(summary(lm10))
print(summary(linmod10))
check.lm(linmod10, lm10, newdata=etit[3:5,])
set.seed(2020)
plotmo(lm10,     pt.col="green", do.par=2)
set.seed(2020)
plotmo(linmod10, pt.col="green", do.par=0)
par(org.par)

cat0("==data.frame with strings\n")

df.with.string <-
    data.frame(1:5,
               c(1,2,-1,4,5),
               c("a", "b", "a", "a", "b"),
               stringsAsFactors=FALSE)
colnames(df.with.string) <- c("num1", "num2", "string")

linmod30.form <- linmod(num1~num2, df.with.string)
lm30       <- lm(num1~num2, df.with.string)
check.lm(linmod30.form, lm30, check.newdata=FALSE)

linmod31.form <- linmod(num1~., df.with.string)
lm31       <- lm(num1~., df.with.string)
check.lm(linmod31.form, lm31, check.newdata=FALSE)

expect.err(try(linmod(string~., df.with.string)), "'y' is not numeric or logical")

vec <- c(1,2,3,4,3)
expect.err(try(linmod(df.with.string, vec)), "non-numeric column in 'x'")
expect.err(try(linmod(etit$pclass, etit$survived)), "non-numeric column in 'x'")

cat0("==x is singular\n")

set.seed(1)
x2 <- matrix(rnorm(6), nrow=2)
y2 <- c(1,2)
expect.err(try(linmod(y2~x2)), "'x' is singular (it has 4 columns but its rank is 2)")

x3 <- matrix(1:10, ncol=2)
y3 <- c(1,2,9,4,5)
expect.err(try(linmod(y3~x3)), "'x' is singular (it has 3 columns but its rank is 2)")

expect.err(try(linmod(trees[1,1:2], trees[1,3])), "'x' is singular (it has 3 columns but its rank is 1)")

x2a <- matrix(1:6, nrow=3)
y2a <- c(1,2,3)
expect.err(try(linmod(y2a~x2a)), "'x' is singular (it has 3 columns but its rank is 2)")

cat0("==perfect fit (residuals are zero)\n")

set.seed(1)
x2b <- matrix(rnorm(6), nrow=3)
y2b <- c(1,2,3)
data.x2b <- data.frame(x2b, y2b)
colnames(data.x2b) <- c("x1", "x2", "y")
linmod.x2b <- linmod(y~., data=data.x2b)
print(summary(linmod.x2b)) # will have "Residual degrees-of-freedom is zero" comment
lm.x2b <- lm(y~., data=data.x2b)
print(summary(lm.x2b)) # will have "ALL 3 residuals are 0" comment
check.lm(linmod.x2b, lm.x2b, newdata=data.x2b[1:2,]+1, check.sigma=FALSE)

x2c <- 1:10
y2c <- 11:20
data.x2c <- data.frame(x2c, y2c)
colnames(data.x2c) <- c("x", "y")
linmod.x2c <- linmod(y~., data=data.x2c)
print(summary(linmod.x2c))
lm.x2c <- lm(y~., data=data.x2c)
options(warn=1) # print warnings as they occur
print(summary(lm.x2c)) # will have "essentially perfect fit: summary may be unreliable" comment
options(warn=2) # treat warnings as errors
check.lm(linmod.x2c, lm.x2c, newdata=data.x2c[1:2,]+1, check.sigma=FALSE)

par(mfrow=c(2,2)) # all plots on same page so can compare
plot(linmod.x2b, main="linmod.x2b\nall residuals are zero")
plot(lm.x2b, which=1, main="lm.x2b")
plot(linmod.x2c, main="linmod.x2c")
plot(lm.x2c, which=1, main="lm.x2c")
par(org.par)

cat0("==nrow(x) does not match length(y)\n")

x4 <- matrix(1:10, ncol=2)
y4 <- c(1,2,9,4)
expect.err(try(linmod(x4, y4)), "nrow(x) is 5 but length(y) is 4")

x5 <- matrix(1:10, ncol=2)
y5 <- c(1,2,9,4,5,9)
expect.err(try(linmod(x5, y5)), "nrow(x) is 5 but length(y) is 6")

cat0("==y has multiple columns\n")

vec <- c(1,2,3,4,3)
y2 <- cbind(c(1,2,3,4,9), vec^2)
expect.err(try(linmod(vec, y2)), "nrow(x) is 5 but length(y) is 10")
expect.err(try(linmod(y2~vec)), "nrow(x) is 5 but length(y) is 10")

cat0("==NA in x\n")

x <- tr[,1:2]
y <- tr[,3]
x[2,2] <- NA
expect.err(try(linmod(x, y)), "NA in 'x'")

x <- tr[,1:2]
y <- tr[,3]
y[9] <- NA
expect.err(try(linmod(x, y)), "NA in 'y'")

# Following added Sep 2020 (prior to this, predict.linmod gave an incorrect error message)
cat0("==test formulas that use functions on rhs variables, like Volume~sqrt(Girth)\n")

linmod.sqrt1 <- linmod(Volume~sqrt(as.numeric(Girth)), data=tr)
cat0("==print(summary(linmod.sqrt1))\n")
print(summary(linmod.sqrt1))
lm.sqrt1 <- lm(Volume~sqrt(as.numeric(Girth)), data=tr)
check.lm(linmod.sqrt1, lm.sqrt1)
stopifnot(almost.equal(predict(linmod.sqrt1, newdata=data.frame(Girth=10, Height=80)),
                       predict(lm.sqrt1,     newdata=data.frame(Girth=10, Height=80))))
stopifnot(almost.equal(predict(linmod.sqrt1, newdata=as.matrix(tr[1:3,])),
                       predict(lm.sqrt1,     newdata=tr[1:3,])))
par(mfrow=c(2,2)) # all plots on same page so can compare
plotmo(linmod.sqrt1, do.par=FALSE)
plotmo(lm.sqrt1,     do.par=FALSE)
par(org.par)

linmod.sqrt2 <- linmod(Volume~sqrt(Girth)+Height+Girth, data=tr)
cat0("==print(summary(linmod.sqrt2))\n")
print(summary(linmod.sqrt2))
lm.sqrt2 <- lm(Volume~sqrt(Girth)+Height+Girth, data=tr)
check.lm(linmod.sqrt2, lm.sqrt2)
stopifnot(almost.equal(predict(linmod.sqrt2, newdata=data.frame(Girth=10, Height=80)),
                       predict(lm.sqrt2,     newdata=data.frame(Girth=10, Height=80))))
stopifnot(almost.equal(predict(linmod.sqrt2, newdata=as.matrix(tr[1:3,])),
                       predict(lm.sqrt2,     newdata=tr[1:3,])))
par(mfrow=c(2,2)) # all plots on same page so can compare
plotmo(linmod.sqrt2, do.par=FALSE)
plotmo(lm.sqrt2,     do.par=FALSE)
par(org.par)

lm.sqrt.as.numeric <- lm(survived~sqrt(age)+as.numeric(pclass), data=etit)
linmod.sqrt.as.numeric <- linmod(survived~sqrt(age)+as.numeric(pclass), data=etit)
check.lm(linmod.sqrt.as.numeric, lm.sqrt.as.numeric, newdata=etit[3:5,])
expect.err(try(predict(linmod.sqrt.as.numeric, newdata=data.frame(age=30))), # pclass missing
           "object 'pclass' not found")
par(mfrow=c(2,2)) # all plots on same page so can compare
plotmo(linmod.sqrt.as.numeric, do.par=FALSE)
plotmo(lm.sqrt.as.numeric,     do.par=FALSE)
par(org.par)

y.age    <- etit[,"age"]
x.pclass <- etit[,"pclass"]
x.sex    <- etit[,"sex"]
linmod.y.age.sex.pclass <- linmod(y.age ~ as.numeric(x.pclass) + x.sex)
lm.y.age.sex.pclass     <- lm(    y.age ~ as.numeric(x.pclass) + x.sex)
stopifnot(identical(linmod.y.age.sex.pclass$coef, lm.y.age.sex.pclass$coef))
options(warn=1) # print warnings as they occur to test stop() in linmod.R::process.newdata.formula
# TODO following says variable 'as.numeric(x.pclass)' may be missing
#      it should say  variable 'x.pclass' may be missing
expect.err(try(predict(linmod.y.age.sex.pclass, newdata=etit[3:5,1,drop=FALSE])),
           "newdata has 3 rows but model.frame returned 18 rows (variable 'as.numeric(x.pclass)' may be missing from newdata)")
options(warn=2) # back to treating warnings as errors

cat0("==misc tests with different kinds of data\n")

data3 <- data.frame(s=c("a", "b", "a", "c", "a"), num=c(1,5,1,9,2), y=c(1,3,2,5,3), stringsAsFactors=F)
stopifnot(sapply(data3, class) == c("character", "numeric", "numeric"))
a40 <- linmod(y~., data=data3)
print(summary(a40))
stopifnot(almost.equal(a40$coefficients, c(0, -4.5, -8.5, 1.5), max=0.001))
stopifnot(almost.equal(predict(a40, newdata=data3[2:3,]),
                       c(3.0, 1.5), max=0.001))

data4 <- data.frame(s=c("a", "b", "a", "c", "a"), num=c(1,5,1,9,2), y=c(1,3,2,5,3), stringsAsFactors=T)
stopifnot(sapply(data4, class) == c("factor", "numeric", "numeric"))
expect.err(try(linmod(data4[,1:2], data4[,3])), "non-numeric column in 'x'")

# following gives no error (and matches lm) even though col 1 of data3 is character not factor
a41 <- linmod(y~., data=data4)
print(summary(a41))
stopifnot(almost.equal(predict(a41, newdata=data3[2:3,]),
                       c(3.0, 1.5), max=0.001))

data5 <- data.frame(s=c("a", "b", "c", "a", "a"), num=c(1,9,4,2,6), y=c(1,2,3,5,3), stringsAsFactors=F)
stopifnot(almost.equal(predict(a41, newdata=data5[1:3,1:2]),
                        c(1.5, 9.0, -2.5), max=0.001))

data6 <- data.frame(s=c("a", "b", "c", "a9", "a"),
                    num=c(1,9,4,2,6),
                    num2=c(1,9,4,2,7),
                    y=c(1,2,3,5,3), stringsAsFactors=T)
expect.err(try(predict(a41, newdata=data6[1:3,1])), "object 's' not found")
expect.err(try(predict(a41, newdata=data6[1:3,c(1,1)])), "object 'num' not found")

expect.err(try(predict(a41, newdata=data.frame(s=1, num=2, y=3))), "variable 's' is not a factor")

expect.err(try(predict(a41, newdata=1:9)),
           "object 's' not found")

expect.err(try(predict(a41, newdata=data.frame())), "'newdata' is empty")

# perfect fit (residuals are all zero)
linmod.data6 <- linmod(y~s+num, data=data6)
print(summary(linmod.data6))
lm.data6 <- lm(y~s+num, data=data6)
print(summary(lm.data6))
check.lm(linmod.data6, lm.data6, newdata=data6[2,,drop=FALSE], check.sigma=FALSE)

expect.err(try(linmod(y~., data=data6)), "'x' is singular (it has 6 columns but its rank is 5)")

tr.na <- trees
tr.na[9,3] <- NA
expect.err(try(linmod(Volume~.,data=tr.na)), "NA in 'y'")
expect.err(try(linmod(tr.na[,1:2], tr.na[,3])), "NA in 'y'")

tr.na <- trees
tr.na[10,1] <- NA
expect.err(try(linmod(Volume~.,data=tr.na)), "NA in 'x'")
expect.err(try(linmod(tr.na[,1:2], tr.na[,3])), "NA in 'x'")

a42 <- linmod(trees[,1:2], trees[, 3])
newdata1 <- data.frame(Girth=20)
expect.err(try(predict(a42, newdata=newdata1)), "ncol(newdata) is 1 but should be 2")
newdata3 <- data.frame(Girth=20, extra1=21, extra2=22)
expect.err(try(predict(a42, newdata=newdata3)), "ncol(newdata) is 3 but should be 2")
expect.err(try(predict(a42, newdata=data.frame())), "'newdata' is empty")
newdata.with.NA <- data.frame(Girth=20, Height=NA)
expect.err(try(predict(a42, newdata=newdata.with.NA)), "NA in 'newdata'")

a43 <- linmod(Volume~.,data=trees)
expect.err(try(predict(a43, newdata=newdata.with.NA)), "NA in 'newdata'")
lm43 <- lm(Volume~.,data=trees)
# message from predict.lm could be better
expect.err(try(predict(lm43, newdata=newdata.with.NA)),
              "variable 'Height' was fitted with type \"numeric\" but type \"logical\" was supplied")

y6 <- 1:5
x6 <- data.frame()
options(warn=1) # print warnings as they occur
expect.err(try(linmod(x6, y6)), "'x' is empty")
options(warn=2) # treat warnings as errors

y7 <- data.frame()
x7 <- 1:5
expect.err(try(linmod(x7, y7)), "'y' is empty")

# duplicated column names
data7 <- matrix(1:25, ncol=5)
colnames(data7) <- c("y", "x1", "x1", "x3", "x4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x1\" in 'x' is duplicated")

colnames(data7) <- c("y", "x1", "x2", "x2", "x4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x2\" in 'x' is duplicated")

colnames(data7) <- c("y", "x1", "x2", "x2", "x2")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"x2\" in 'x' is duplicated")

# column name V2 will be created but it clashes with the existing column name
colnames(data7) <- c("y", "V2", "", "V3", "V4")
expect.err(try(linmod(data7[,-1], data7[,1])), "column name \"V2\" in 'x' is duplicated")

# missing column names
trees1 <- trees
colnames(trees1) <- NULL
cat0("a52\n")
a52 <- linmod(trees1[,1:2], trees1[,3])
print(summary(a52))

trees1 <- trees
colnames(trees1) <- c("", "Height", "Volume") # was Girth Height Volume
cat0("linmod.form.Volume.trees1\n")
linmod.form.Volume.trees1 <- linmod(trees1[,1:2], trees1[,3])
print(summary(linmod.form.Volume.trees1))
cat0("linmod.form.Volume.trees1.formula\n")
expect.err(try(linmod(Volume~., data=trees1)), "attempt to use zero-length variable name")

# very long names to test formatting in summary.linmod
trees1 <- trees
colnames(trees1) <- c("Girth.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name",
                      "Height.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name",
                      "Volume.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name")
cat0("a55\n")
a55 <- linmod(Volume.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name~
              Girth.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name+
              Height.a.very.long.name.in.fact.an.exceptionally.exceptionally.exceptionally.long.name,
              data=trees1)
print(summary(a55))

# intercept-only model
intonly.form <- linmod(Volume~1, data=trees)
print(summary(intonly.form))
stopifnot(length(coef(intonly.form)) == 1)
try(plotmo(intonly.form)) # Error in plotmo(intonly.form) : x is empty
plotres(intonly.form)
expect.err(try(plotmo(intonly.form)), "x is empty")
expect.err(try(linmod(rep(1, length.out=nrow(trees)), trees$Volume)),
              "'x' is singular (it has 2 columns but its rank is 1)")

# various tests for bad args
expect.err(try(linmod(trees[,1:2])), "no 'y' argument")

# test stop.if.dot.arg.used
expect.err(try(linmod(Volume~., data=trees, nonesuch=99)),
              "unused argument (nonesuch = 99)")
expect.err(try(linmod(trees[,1:2], trees[,3], nonesuch=linmod)),
              "unused argument (nonesuch = function (...)")
expect.err(try(summary(linmod(trees[,1:2], trees[,3]), nonesuch=linmod)),
              "unused argument (nonesuch = function (...)")
expect.err(try(print(linmod(trees[,1:2], trees[,3]), nonesuch=linmod)),
              "unused argument (nonesuch = function (...)")
expect.err(try(predict(linmod.form.Volume.tr, nonesuch=99)),
              "unused argument (nonesuch = 99)")

# check partial matching on type argument
stopifnot(identical(predict(linmod.form.Volume.tr, type="r"),    predict(linmod.form.Volume.tr)))
stopifnot(identical(predict(linmod.form.Volume.tr, type="resp"), predict(linmod.form.Volume.tr)))
expect.err(try(predict(linmod.form.Volume.tr, type="nonesuch")), "'arg' should be \"response\"")

# test additional method functions (see linmod.methods.R)

check.lm(linmod.form.Volume.tr, lm.Volume.tr, newdata=trees[3,1:2])
stopifnot(almost.equal(coef(linmod.form.Volume.tr), coef(lm.Volume.tr)))
stopifnot(identical(names(coef(linmod.form.Volume.tr)), names(coef(lm.Volume.tr))))
stopifnot(almost.equal(fitted(linmod.form.Volume.tr), fitted(lm.Volume.tr)))
stopifnot(identical(names(fitted(linmod.form.Volume.tr)), names(fitted(lm.Volume.tr))))
stopifnot(identical(na.action(linmod.form.Volume.tr), na.action(lm.Volume.tr)))
stopifnot(almost.equal(residuals(linmod.form.Volume.tr), residuals(lm.Volume.tr)))
stopifnot(identical(names(residuals(linmod.form.Volume.tr)), names(residuals(lm.Volume.tr))))
stopifnot(identical(names(case.names(linmod.form.Volume.tr)), names(case.names(lm.Volume.tr))))
stopifnot(identical(variable.names(linmod.form.Volume.tr), variable.names(lm.Volume.tr)))
stopifnot(identical(nobs(linmod.form.Volume.tr), nobs(lm.Volume.tr)))
stopifnot(identical(weights(linmod.form.Volume.tr), weights(lm.Volume.tr)))
stopifnot(almost.equal(df.residual(linmod.form.Volume.tr), df.residual(lm.Volume.tr)))
stopifnot(identical(names(df.residual(linmod.form.Volume.tr)), names(df.residual(lm.Volume.tr))))
stopifnot(almost.equal(deviance(linmod.form.Volume.tr), deviance(lm.Volume.tr)))
stopifnot(identical(names(deviance(linmod.form.Volume.tr)), names(deviance(lm.Volume.tr))))
stopifnot(identical(weights(linmod.form.Volume.tr), weights(lm.Volume.tr)))
stopifnot(identical(model.frame(linmod.form.Volume.tr), model.frame(lm.Volume.tr)))
stopifnot(identical(model.matrix(linmod.form.Volume.tr), model.matrix(lm.Volume.tr)))
stopifnot(identical(model.matrix(linmod.form.Volume.tr, data=tr[1:2,]),
                    model.matrix(lm.Volume.tr,          data=tr[1:2,])))
stopifnot(almost.equal(logLik(linmod.form.Volume.tr), logLik(lm.Volume.tr)))
expect.err(try(logLik(linmod.form.Volume.tr, REML=TRUE)), "!REML is not TRUE")
library(sandwich) # for estfun.lm
stopifnot(almost.equal(estfun(linmod.form.Volume.tr), estfun(lm.Volume.tr)))

linmod.form.Volume.tr.update <- update(linmod.form.Volume.tr, formula.=Volume~Height)
lm.Volume.tr.update          <- update(lm.Volume.tr, formula.=Volume~Height)
check.lm(linmod.form.Volume.tr.update, lm.Volume.tr.update)

check.lm(linmod.xy.Volume.tr, lm.Volume.tr, newdata=trees[3,1:2])
stopifnot(almost.equal(coef(linmod.xy.Volume.tr), coef(lm.Volume.tr)))
stopifnot(identical(names(coef(linmod.xy.Volume.tr)), names(coef(lm.Volume.tr))))
stopifnot(almost.equal(fitted(linmod.xy.Volume.tr), fitted(lm.Volume.tr)))
stopifnot(identical(names(fitted(linmod.xy.Volume.tr)), names(fitted(lm.Volume.tr))))
stopifnot(identical(na.action(linmod.xy.Volume.tr), na.action(lm.Volume.tr)))
stopifnot(almost.equal(residuals(linmod.xy.Volume.tr), residuals(lm.Volume.tr)))
stopifnot(identical(names(residuals(linmod.xy.Volume.tr)), names(residuals(lm.Volume.tr))))
stopifnot(identical(case.names(linmod.xy.Volume.tr), case.names(lm.Volume.tr)))
stopifnot(identical(variable.names(linmod.xy.Volume.tr), variable.names(lm.Volume.tr)))
stopifnot(identical(nobs(linmod.xy.Volume.tr), nobs(lm.Volume.tr)))
stopifnot(identical(weights(linmod.xy.Volume.tr), weights(lm.Volume.tr)))
stopifnot(almost.equal(df.residual(linmod.xy.Volume.tr), df.residual(lm.Volume.tr)))
stopifnot(identical(names(df.residual(linmod.xy.Volume.tr)), names(df.residual(lm.Volume.tr))))
stopifnot(almost.equal(deviance(linmod.xy.Volume.tr), deviance(lm.Volume.tr)))
stopifnot(identical(names(deviance(linmod.xy.Volume.tr)), names(deviance(lm.Volume.tr))))
stopifnot(identical(weights(linmod.xy.Volume.tr), weights(lm.Volume.tr)))
expect.err(try(model.frame(linmod.xy.Volume.tr)),  "model.frame cannot be used on linmod models built without a formula")
expect.err(try(model.matrix(linmod.xy.Volume.tr)),
              "model.frame cannot be used on linmod models built without a formula")
stopifnot(almost.equal(logLik(linmod.xy.Volume.tr), logLik(lm.Volume.tr)))

par(mfrow=c(2,2))
plot(linmod.form.Volume.tr)
plot(lm.Volume.tr, which=1, main="lm.Volume.tr")
plot(linmod.xy.Volume.tr)
plot(linmod.form.Volume.tr, xlim=c(0,80), ylim=c(-10,10), pch=20, main="linmod.form.Volume.tr: test plot args")
par(org.par)

cat0("==test one predictor model\n")

linmod.onepred.form <- linmod(Volume~Girth, data=tr) # one predictor
lm.onepred.form <- lm(Volume~Girth, data=tr)
check.lm(linmod.onepred.form, lm.onepred.form, newdata=trees[3,1:2])
linmod.onepred.xy <- linmod(tr[,1,drop=FALSE], tr[,3]) # one predictor
print(summary(linmod.onepred.xy))
check.lm(linmod.onepred.xy, lm.onepred.form, newdata=trees[3,1,drop=FALSE])

par(mfrow=c(2,2))
plot(linmod.onepred.form)
plot(lm.onepred.form, which=1, main="lm.onepred.form")
plot(linmod.onepred.xy)
par(org.par)
plotres(linmod.onepred.form)
plotmo(linmod.onepred.form, pt.col=2)

cat0("==test no intercept model\n")
# no intercept models are only supported with the formula interface (not x,y interface)

linmod.noint <- linmod(Volume~.-1, data=trees) # no intercept
print(summary(linmod.noint))
lm.noint <- lm(Volume~.-1, data=trees) # no intercept
check.lm(linmod.noint, lm.noint)
linmod.noint.keep <- linmod(Volume~.-1, data=trees, keep=TRUE)
print(summary(linmod.noint.keep))

check.lm(linmod.noint, lm.noint)
stopifnot(class(linmod.noint.keep$data)   == class(linmod.form.Volume.trees.keep$data))
stopifnot(all(dim(linmod.noint.keep$data) == dim(linmod.form.Volume.trees.keep$data)))
stopifnot(all(linmod.noint.keep$data == linmod.form.Volume.trees.keep$data))
stopifnot(class(linmod.noint.keep$y)   == class(linmod.form.Volume.trees.keep$y))
stopifnot(all(dim(linmod.noint.keep$data) == dim(linmod.form.Volume.trees.keep$data)))
stopifnot(all(linmod.noint.keep$data == linmod.form.Volume.trees.keep$data))

# check method functions in no-intercept model
stopifnot(almost.equal(coef(linmod.noint), coef(lm.noint)))
stopifnot(identical(names(coef(linmod.noint)), names(coef(lm.noint))))
stopifnot(almost.equal(fitted(linmod.noint), fitted(lm.noint)))
stopifnot(identical(names(fitted(linmod.noint)), names(fitted(lm.noint))))
stopifnot(identical(na.action(linmod.noint), na.action(lm.noint)))
stopifnot(almost.equal(residuals(linmod.noint), residuals(lm.noint)))
stopifnot(identical(names(residuals(linmod.noint)), names(residuals(lm.noint))))
stopifnot(identical(case.names(linmod.noint), case.names(lm.noint)))
stopifnot(identical(variable.names(linmod.noint), variable.names(lm.noint)))
stopifnot(identical(nobs(linmod.noint), nobs(lm.noint)))
stopifnot(identical(weights(linmod.noint), weights(lm.noint)))
stopifnot(almost.equal(df.residual(linmod.noint), df.residual(lm.noint)))
stopifnot(identical(names(df.residual(linmod.noint)), names(df.residual(lm.noint))))
stopifnot(almost.equal(deviance(linmod.noint), deviance(lm.noint)))
stopifnot(identical(names(deviance(linmod.noint)), names(deviance(lm.noint))))
stopifnot(identical(weights(linmod.noint), weights(lm.noint)))
stopifnot(identical(model.frame(linmod.noint), model.frame(lm.noint)))
stopifnot(identical(model.matrix(linmod.noint), model.matrix(lm.noint)))
stopifnot(identical(model.matrix(linmod.noint, data=tr[1:2,]),
                    model.matrix(lm.noint,     data=tr[1:2,])))
stopifnot(almost.equal(logLik(linmod.noint), logLik(lm.noint)))
stopifnot(almost.equal(estfun(linmod.noint), estfun(lm.noint)))

# check error messages with bad newdata in no-intercept model
expect.err(try(predict(linmod.noint, newdata=NA)),
              "object 'Girth' not found")
expect.err(try(predict(linmod.noint, newdata=data.frame(Height=c(1,NA), Girth=c(3,4)))),
              "NA in 'newdata'")
expect.err(try(predict(linmod.noint, newdata=trees[0,])), "'newdata' is empty")
expect.err(try(predict(linmod.noint, newdata=trees[3:5,"Height"])), "object 'Girth' not found")
# check that extra fields in predict newdata are ok with (formula) models without intercept
stopifnot(almost.equal(predict(linmod.noint, newdata=data.frame(Girth=10, Height=80, extra=99)),
                       predict(lm.noint,  newdata=data.frame(Girth=10, Height=80, extra=99))))

par(mfrow=c(2,2))
plot(linmod.noint)
plot(lm.noint, which=1, main="lm.noint")
par(org.par)

plotres(linmod.noint)
plotmo(linmod.noint)

cat0("==test one predictor no intercept model\n")
# no intercept models are only supported with the formula interface (not x,y interface)

linmod.onepred.noint <- linmod(Volume~Girth-1, data=trees) # one predictor, no intercept
print(summary(linmod.onepred.noint))
lm.onepred.noint <- lm(Volume~Girth-1, data=trees) # one predictor, no intercept
check.lm(linmod.onepred.noint, lm.onepred.noint)
linmod.onepred.noint.keep <- linmod(Volume~.-1, data=trees, keep=TRUE)
print(summary(linmod.onepred.noint.keep))

check.lm(linmod.onepred.noint, lm.onepred.noint)
stopifnot(class(linmod.onepred.noint.keep$data)   == class(linmod.form.Volume.trees.keep$data))
stopifnot(all(dim(linmod.onepred.noint.keep$data) == dim(linmod.form.Volume.trees.keep$data)))
stopifnot(all(linmod.onepred.noint.keep$data == linmod.form.Volume.trees.keep$data))
stopifnot(class(linmod.onepred.noint.keep$y)   == class(linmod.form.Volume.trees.keep$y))
stopifnot(all(dim(linmod.onepred.noint.keep$data) == dim(linmod.form.Volume.trees.keep$data)))
stopifnot(all(linmod.onepred.noint.keep$data == linmod.form.Volume.trees.keep$data))

# check method functions in one predictor no-intercept model
stopifnot(almost.equal(coef(linmod.onepred.noint), coef(lm.onepred.noint)))
stopifnot(identical(names(coef(linmod.onepred.noint)), names(coef(lm.onepred.noint))))
stopifnot(almost.equal(fitted(linmod.onepred.noint), fitted(lm.onepred.noint)))
stopifnot(identical(names(fitted(linmod.onepred.noint)), names(fitted(lm.onepred.noint))))
stopifnot(identical(na.action(linmod.onepred.noint), na.action(lm.onepred.noint)))
stopifnot(almost.equal(residuals(linmod.onepred.noint), residuals(lm.onepred.noint)))
stopifnot(identical(names(residuals(linmod.onepred.noint)), names(residuals(lm.onepred.noint))))
stopifnot(identical(case.names(linmod.onepred.noint), case.names(lm.onepred.noint)))
stopifnot(identical(variable.names(linmod.onepred.noint), variable.names(lm.onepred.noint)))
stopifnot(identical(nobs(linmod.onepred.noint), nobs(lm.onepred.noint)))
stopifnot(identical(weights(linmod.onepred.noint), weights(lm.onepred.noint)))
stopifnot(almost.equal(df.residual(linmod.onepred.noint), df.residual(lm.onepred.noint)))
stopifnot(identical(names(df.residual(linmod.onepred.noint)), names(df.residual(lm.onepred.noint))))
stopifnot(almost.equal(deviance(linmod.onepred.noint), deviance(lm.onepred.noint)))
stopifnot(identical(names(deviance(linmod.onepred.noint)), names(deviance(lm.onepred.noint))))
stopifnot(identical(weights(linmod.onepred.noint), weights(lm.onepred.noint)))
stopifnot(identical(model.frame(linmod.onepred.noint), model.frame(lm.onepred.noint)))
stopifnot(identical(model.matrix(linmod.onepred.noint), model.matrix(lm.onepred.noint)))
stopifnot(identical(model.matrix(linmod.onepred.noint, data=tr[1:2,]),
                    model.matrix(lm.onepred.noint,     data=tr[1:2,])))
stopifnot(almost.equal(logLik(linmod.onepred.noint), logLik(lm.onepred.noint)))
stopifnot(almost.equal(estfun(linmod.onepred.noint), estfun(lm.onepred.noint)))

# check error messages with bad newdata in one predictor no-intercept model
expect.err(try(predict(linmod.onepred.noint, newdata=99)), "object 'Girth' not found")
expect.err(try(predict(linmod.onepred.noint, newdata=data.frame(Girth=NA))), "NA in 'newdata'")
expect.err(try(predict(linmod.onepred.noint, newdata=trees[0,1])), "'newdata' is empty")
expect.err(try(predict(linmod.onepred.noint, newdata=trees[3:5,"Height"])), "object 'Girth' not found")
# check that extra fields in predict newdata are ok with (formula) models without intercept
stopifnot(almost.equal(predict(linmod.onepred.noint, newdata=data.frame(Girth=10, extra=99)),
                       predict(lm.onepred.noint,     newdata=data.frame(Girth=10, extra=99))))

par(mfrow=c(2,2))
plot(linmod.onepred.noint)
plot(lm.onepred.noint, which=1, main="lm.onepred.noint")
par(org.par)

plotres(linmod.onepred.noint)
plotmo(linmod.onepred.noint)

expect.err(try(linmod(Volume~nonesuch, data=trees)), "object 'nonesuch' not found")
expect.err(try(linmod(Volume~0, data=trees)),   "'x' is empty") # no predictor
expect.err(try(linmod(Volume~-1, data=trees)), "'x' is empty") # no predictor, no intercept

cat0("==check model with many variables\n")

set.seed(2018)
p <- 300 # number of variables
n <- floor(1.1 * p)
bigdat <- as.data.frame(matrix(rnorm(n * (p+1)), ncol=p+1))
colnames(bigdat) <- c("y", paste0("var", 1:p))
lm.bigdat <- lm(y~., data=bigdat)
linmod.bigdat <- linmod(y~., data=bigdat)
check.lm(linmod.form.Volume.tr, lm.Volume.tr)
print(linmod.bigdat)
print(summary(linmod.bigdat))
expect.err(try(predict(linmod.bigdat, newdata=bigdat[,1:(p-3)])), "object 'var297' not found")
plot(linmod.bigdat)
# plotmo(linmod.bigdat) # works, but commented out because slow(ish)
# plotres(linmod.bigdat) # ditto

cat0("==check use of matrix as data in linmod.form\n")
# linmod.form allows a matrix, lm doesn't TODO is this inconsistency what we want?
tr.mat <- as.matrix(tr)
cat0("class(tr.mat)=", class(tr.mat), "\n") # class(tr.mat)=matrix
expect.err(try(lm(Volume~., data=tr.mat)), "'data' must be a data.frame, not a matrix or an array")
linmod.form.Volume.mat.tr <- linmod(Volume~., data=tr.mat)
check.lm(linmod.form.Volume.mat.tr, linmod.form.Volume.tr)
cat0("==print(summary(linmod.form.Volume.mat.tr))\n")
print(summary(linmod.form.Volume.mat.tr))
plotres(linmod.form.Volume.mat.tr)

tr.mat.no.colnames <- as.matrix(tr)
colnames(tr.mat.no.colnames) <- NULL
expect.err(try(linmod(Volume~., data=tr.mat.no.colnames)), "object 'Volume' not found")
linmod.form.Volume.mat.tr.no.colnames <- linmod(V3~., data=tr.mat.no.colnames)
check.lm(linmod.form.Volume.mat.tr.no.colnames, linmod.form.Volume.tr,
         check.coef.names=FALSE, check.newdata=FALSE) # no check.newdata else object 'V1' not found

# Check what happens when we change the original data used to build the model.
# Use plotres as an example function that must figure out residuals from predict().

pr <- function(model, main=deparse(substitute(model)))
{
    plotres(model, which=3, main=main) # which=3 for just the residuals plot
}
cat0("==linmod.formula: change data used to build the model\n")

trees1 <- trees
linmod.trees1 <- linmod(Volume~., data=trees1)
# delete the saved residuals and fitted.values so plotres has to use the saved
# call etc. to get the x and y used to build the model, and rely on predict()
linmod.trees1$residuals <- NULL
linmod.trees1$fitted.values <- NULL
par(mfrow=c(3,3))
pr(linmod.trees1)
trees1 <- trees[, 3:1]                      # change column order in original data
pr(linmod.trees1, "change col order")
trees1 <- trees[1:3, ]                      # change number of rows in original data
pr(linmod.trees1, "change nbr rows")        # TODO wrong residuals! (lm has the same issue)
cat("call$data now refers to the changed data:\n") # lm has the same problem if called with model=FALSE
print(eval(linmod.trees1$call$data))
cat("model.frame now returns the changed data:\n")
print(model.frame(linmod.trees1))
trees1 <- trees[nrow(tr):1, ]               # change row order (but keep same nbr of rows)
pr(linmod.trees1, "change row order")
colnames(trees1) <- c("x1", "x2", "x3")     # change column names in original data
expect.err(try(pr(linmod.trees1,
              "change colnames")), "cannot get the original model predictors")
trees1 <- "garbage"
expect.err(try(pr(linmod.trees1,
              "trees1=\"garbage\"")), "cannot get the original model predictors")
trees1 <- 1:1000
expect.err(try(pr(linmod.trees1,
              "trees1=1:1000")), "cannot get the original model predictors")
trees1 <- NULL                              # original data no longer available
expect.err(try(pr(linmod.trees1,
              "trees1=NULL")), "cannot get the original model predictors")
remove(trees1)
expect.err(try(pr(linmod.trees1,
              "remove(trees1)")), "cannot get the original model predictors")

# similar to above, but don't delete the saved residuals and fitted.values
trees1 <- trees
linmod2.trees1 <- linmod(Volume~., data=trees1)
trees1 <- trees[1:3, ]                      # change number of rows in original data
expect.err(try(plotmo(linmod2.trees1)), "plotmo_y returned the wrong length (got 3 but expected 31)")

par(org.par)

cat0("==linmod.formula(keep=TRUE): change data used to build the model\n")
par(mfrow=c(3,3))
trees1 <- trees
linmod.trees1.keep <- linmod(Volume~., data=trees1, keep=TRUE)
# delete the saved residuals and fitted.values so plotres has to use the saved
# call etc. to get the x and y used to build the model, and rely on predict()
linmod.trees1.keep$residuals <- NULL
linmod.trees1.keep$fitted.values <- NULL
pr(linmod.trees1.keep)
trees1 <- trees[, 3:1]                      # change column order in original data
pr(linmod.trees1.keep, "change col order")
trees1 <- trees[1:3, ]                      # change number of rows in original data
pr(linmod.trees1.keep, "change nbr rows")
trees1 <- trees[nrow(tr):1, ]               # change row order (but keep same nbr of rows)
pr(linmod.trees1.keep, "change row order")
colnames(trees1) <- c("x1", "x2", "x3")     # change column names in original data
pr(linmod.trees1.keep, "change colnames")
trees1 <- NULL                              # original data no longer available
pr(linmod.trees1.keep, "trees1=NULL")
remove(trees1)
pr(linmod.trees1.keep, "remove(trees1)")
par(org.par)

cat0("==linmod.default: change data used to build the model\n")
trees1 <- trees
x1 <- trees1[,1:2]
y1 <- trees1[,3]
linmod.xy <- linmod(x1, y1)
# delete the saved residuals and fitted.values so plotres has to use the saved
# call etc. to get the x1 and y1 used to build the model, and rely on predict()
linmod.xy$residuals <- NULL
linmod.xy$fitted.values <- NULL
par(mfrow=c(3,3))
pr(linmod.xy)
x1 <- trees1[,2:1]                 # change column order in original x1
pr(linmod.xy, "change col order")
x1 <- trees1[1:3, 1:2]                      # change number of rows in original x1
expect.err(try(pr(linmod.xy, "change nbr rows")),
              "plotmo_y returned the wrong length (got 31 but expected 3)") # TODO different behaviour to linmod.trees1
cat("call$x1 now refers to the changed x1:\n") # lm has the same problem if called with model=FALSE
print(eval(linmod.xy$call$x1))
x1 <- trees1[nrow(tr):1, 1:2]               # change row order (but keep same nbr of rows)
pr(linmod.xy, "change row order")
x1 <- trees1[,1:2]
colnames(x1) <- c("x1", "x2")     # change column names in original x1
pr(linmod.xy, "change colnames")
x1 <- "garbage"
expect.err(try(pr(linmod.xy, "x1=\"garbage\"")), "cannot get the original model predictors")
x1 <- 1:1000
expect.err(try(pr(linmod.xy, "x1=1:1000")), "ncol(newdata) is 1 but should be 2")
x1 <- NULL                              # original x1 no longer available
expect.err(try(pr(linmod.xy, "x1=NULL")), "cannot get the original model predictors")
remove(x1)
expect.err(try(pr(linmod.xy, "remove(x1)")), "cannot get the original model predictors")

# similar to above, but don't delete the saved residuals and fitted.values
trees1 <- trees
x1 <- trees1[,1:2]
y1 <- trees1[,3]
linmod.xy <- linmod(x1, y1)
x1 <- trees1[1:3, 1:2]                      # change number of rows in original x1
expect.err(try(plotmo(linmod2.x1)), "object 'linmod2.x1' not found") # TODO error message misleading?

par(org.par)

cat0("==linmod.default(keep=TRUE): change data used to build the model\n")
par(mfrow=c(3,3))
trees1 <- trees
x1 <- trees1[,1:2]
linmod.xy <- linmod(x1, y1, keep=TRUE)
# delete the saved residuals and fitted.values so plotres has to use the saved
# call etc. to get the x1 and y1 used to build the model, and rely on predict()
linmod.xy$residuals <- NULL
linmod.xy$fitted.values <- NULL
pr(linmod.xy.keep)
x1 <- trees1[, 2:1]                 # change column order in original x1
pr(linmod.xy.keep, "change col order")
x1 <- trees1[1:3, 1:2]                      # change number of rows in original x1
pr(linmod.xy.keep, "change nbr rows")
x1 <- trees1[nrow(tr):1, 1:2]               # change row order (but keep same nbr of rows)
pr(linmod.xy.keep, "change row order")
x1 <- trees1[,1:2]
colnames(x1) <- c("x1", "x2")     # change column names in original x1
pr(linmod.xy.keep, "change colnames")
x1 <- NULL                              # original x1 no longer available
pr(linmod.xy.keep, "x1=NULL")
remove(x1)
pr(linmod.xy.keep, "remove(x1)")
par(org.par)

cat("==test processing a model created in a function with local data\n")

# pr <- function(model, main=deparse(substitute(model)))
# {
#     plotmo(model, degree1=1, degree2=0, pt.col=2, do.par=FALSE, main=main)
# }
pr <- function(model, main=deparse(substitute(model)))
{
    plotres(model, which=3, main=main) # which=3 for just the residuals plot
}
lm.form.func <- function(keep=FALSE)
{
    local.tr <- trees[1:20,]
    lm(Volume~., data=local.tr, model=keep)
}
linmod.form.func <- function(keep=FALSE)
{
    local.tr <- trees[1:20,]
    model <- linmod(Volume~., data=local.tr, keep=keep)
    # delete the saved residuals and fitted.values so plotres has to use the saved
    # call etc. to get the x and y used to build the model, and rely on predict()
    model$residuals <- NULL
    model$fitted.values <- NULL
    model
}
linmod.xy.func <- function(keep)
{
    xx <- trees[1:20,1:2]
    yy <- trees[1:20,3]
    model <- linmod(xx, yy, keep=keep)
    # delete the saved residuals and fitted.values so plotres has to use the saved
    # call etc. to get the x and y used to build the model, and rely on predict()
    model$residuals <- NULL
    model$fitted.values <- NULL
    model
}
par(mfrow=c(3,2))

lm.form <- lm.form.func(keep=FALSE)
pr(lm.form)

lm.form.keep <- lm.form.func(keep=TRUE)
pr(lm.form.keep)

linmod.form <- linmod.form.func(keep=FALSE)
pr(linmod.form)

linmod.form.keep <- linmod.form.func(keep=TRUE)
pr(linmod.form.keep)

linmod.xy <- linmod.xy.func(keep=FALSE)
expect.err(try(pr(linmod.xy)), "cannot get the original model predictors")

linmod.xy.keep <- linmod.xy.func(keep=TRUE)
pr(linmod.xy.keep)

par(org.par)

# test xlevels (predict with newdata using a string to represent a factor)
data(iris)
linmod.Sepal.Length <- linmod(Sepal.Length~Species,data=iris)
lm.Sepal.Length     <- lm(Sepal.Length~Species,data=iris)
predict.linmod <- predict(linmod.Sepal.Length, newdata=data.frame(Species="setosa"))
predict.lm     <- predict(lm.Sepal.Length,     newdata=data.frame(Species="setosa"))
stopifnot(all.equal(predict.linmod, predict.lm))

source("test.epilog.R")
