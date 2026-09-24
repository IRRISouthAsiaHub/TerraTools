# test.unusual.vars.R: test unusual variable names, and unusual formulas
#
# This file was initially created for plotmo 3.6.0 (Sep 2020)
# ALso tests the naken() func introduced in plotmo 3.6.0 and earth 5.2.0 (Sep 2020)

source("test.prolog.R")
library(earth)
data(ozone1)
data(etitanic)
options(warn=1) # print warnings as they occur

check.naken <- function(s, expected, trace=0)
{
    nude <- plotmo:::naken.formula.string(s, trace=trace)
    printf("%-60.60s %-s\n", s, nude)
    stopifnot(nude == expected)
}
printf("=== check naken.formula.string\n")

# edge cases
check.naken("", "")
check.naken(" ", "")
check.naken("y~", "y ~ ")
check.naken("y~ ", "y ~ ")
check.naken("y ~ ", "y ~ ")
check.naken("y ~  ", "y ~ ")
check.naken(" y ~  ", "y ~ ")
check.naken("[", "[", trace=2)
check.naken("`", "`", trace=2)
# standard formulas
check.naken("x", "x")
check.naken("x1", "x1")
check.naken("y ~ x1 : x2 + x3", "y ~ x1 + x2 + x3", trace=2)
check.naken("y ~ x1 + x2 - x3", "y ~ x1 + x2 + x3", trace=2) # TODO "-" is treated as a "+"
check.naken("y ~ .-x3", "y ~ . + x3")
check.naken("cbind(damage, 6-damage)~temp", "cbind(damage, 6-damage) ~ temp", trace=2)
check.naken("depIndex~q_4 + q_2102+q_2104  +q_3105+  q_3106", "depIndex ~ q_4 + q_2102 + q_2104 + q_3105 + q_3106")
check.naken("doy ~ (vh+wind+humidity)^2", "doy ~ vh + wind + humidity")
check.naken("doy ~ s(wind) + s(humidity,wind) + s(vh)", "doy ~ wind + humidity + vh")
check.naken("log(doy) ~ I(vh*wind) + I(humidity*temp)+log(doy)", "log(doy) ~ vh + wind + humidity + temp + doy")
check.naken("log(doy)~vh+wind+humidity+I(wind*humidity)+temp+log(ibh)", "log(doy) ~ vh + wind + humidity + temp + ibh", trace=2)
check.naken("O3 ~ s(humidity)+s(temp)+s(ibt)+s(temp,ibt)", "O3 ~ humidity + temp + ibt")
check.naken("Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp)", "Ozone^(1/3) ~ Solar.R + Wind + Temp")
check.naken("Volume~(Girth*Height2)-Height", "Volume ~ Girth + Height2 + Height")
check.naken("y ~ s(x) + s(x,z1)", "y ~ x + z1")
check.naken("y~s(x0,x1,k=12)+s(x2)+s(x3,k=20,fx=20)", "y ~ x0 + x1 + x2 + x3")
check.naken("y~x[,1]+x[,2]", "y ~ x[,1] + x[,2]")
check.naken("y~x[,1]+x[,my.list$j]", "y ~ x[,1] + x[,my.list$j]")
check.naken("y~x[,i]+x[,2]", "y ~ x[,i] + x[,2]")
check.naken("Salary~Hitters[,1]", "Salary ~ Hitters[,1]", trace=2)
check.naken("Salary~Hitters[,-1]", "Salary ~ Hitters[,-1]", trace=2)
check.naken("Salary~Hitters[,c(1,2)]", "Salary ~ Hitters[,c(1,2)]", trace=2)
check.naken("Salary~Hitters[,1:2]", "Salary ~ Hitters[,1:2]")
check.naken("Salary~Hitters[,c(1,2)]", "Salary ~ Hitters[,c(1,2)]", trace=2)
# nested brackets
check.naken("y ~ x1[[2]] + x1[[3]]", "y ~ x1[[2]] + x1[[3]]")
check.naken("y[ , 1 ] ~ x1[[2]]", "y[ , 1 ] ~ x1[[2]]")
check.naken("y ~ x0[,nonesuch1 + x1[nsuch2^2  +   3  ]]", "y ~ x0[,nonesuch1 + x1[nsuch2^2  +   3  ]]")
check.naken("y ~ x0[1,x2[3]] + x4[[5]] + x6[ x7[, 8], x9[ ,x10[11] ], drop=x12[13]]", "y ~ x0[1,x2[3]] + x4[[5]] + x6[ x7[, 8], x9[ ,x10[11] ], drop=x12[13]]")
# backquotes
check.naken("y ~ `a b c10` + `def`", "y ~ `a b c10` + `def`")
check.naken("`y` ~ `a b c10` + `def` + s(sqrt(`x 1`))", "`y` ~ `a b c10` + `def` + `x 1`")
# without a response
check.naken("x1 + x[,1] + `x3`", "x1 + x[,1] + `x3`")
check.naken("Salary~Hitters[,c(1,2)]+sqrt(x)", "Salary ~ Hitters[,c(1,2)] + x")
check.naken("Salary~Hitters[,c(1,2)]+sqrt(x)+x99", "Salary ~ Hitters[,c(1,2)] + x + x99")
check.naken("Salary~x1+x2+`x6`+x3", "Salary ~ x1 + x2 + `x6` + x3")
check.naken("x[,c(1,2)] + x[,3]", "x[,c(1,2)] + x[,3]")
check.naken("x[,1] + x[,2] + x[,3] + x[,29] + x[,-14]", "x[,1] + x[,2] + x[,3] + x[,29] + x[,-14]")
check.naken("x[,c(1,2)] + x[,3] + x[,5:6] + x[,-1]", "x[,c(1,2)] + x[,3] + x[,5:6] + x[,-1]")
check.naken("log(y) ~ x9 + ns(x2,4) + s(x3,x4,df=4) + x5:sqrt(x6)", "log(y) ~ x9 + x2 + x3 + x4 + x5 + x6")
check.naken("log(y) ~ x9 + sqrt(x6) + ns(x2,4) + s(x3,x4,df=4) + x5", "log(y) ~ x9 + x6 + x2 + x3 + x4 + x5")
check.naken("x[,1] + sqrt(x2) + 2.34e6 + 1", "x[,1] + x2 + 1")

printf("\n=== test problem in lm() formula with -nonesuch ===\n")

# Using "-nonesuch" in a "." formula (where nonesuch is a non-existent variable name)
# causes the following error in stats::terms.formula (called via model.frame.default)
# Error in terms.formula(formula, data = data) :  (converted from warning)
#   'varlist' has changed (from nvar=3) to new 4 after EncodeVars() -- should no longer happen!
options(warn=2) # treat warnings as errors
expect.err(try(lm(formula = Volume ~ . - nonesuch, data=trees)),
           "'varlist' has changed (from nvar=3) to new 4 after EncodeVars() -- should no longer happen!")
options(warn=1) # print warnings as they occur

printf("\n=== test variables names with spaces in them ===\n")
spaced.trees <- trees
stopifnot(colnames(spaced.trees) == c("Girth", "Height", "Volume")) # sanity check
colnames(spaced.trees) <- c("Girth extra", "Height 999", "Volume") # put spaces in the names

lm.spaced.trees <- lm(Volume~., data=spaced.trees)
options(warn=2)
expect.err(try(plotmo(lm.spaced.trees)),
           "Cannot determine which variables to plot in degree2 plots")
options(warn=1)
plotmo(lm.spaced.trees)            # warning, but still plots (no degree2 plots)
plotmo(lm.spaced.trees, all2=TRUE) # no warning

earth.spaced.trees <- earth(Volume~. , data=spaced.trees, degree=2)
plotmo(earth.spaced.trees)
cat("\nevimp(earth.spaced.trees)\n")
print(evimp(earth.spaced.trees))

printf("\n=== test non standard variable names and use of earth's bx matrix ===\n")
emod <- earth(survived~., data=etitanic, degree=2)
plotmo(emod)
cat("\nevimp(emod)\n")
print(evimp(emod))
bx <- emod$bx
bx.df <- as.data.frame(bx[,-1])  # -1 to drop intercept
bx.df$survived <- etitanic$survived
# following gsub make it a bit easier to see what's going on
# because the next call to earth also creates hinge functions
# (so we end up with nested hinge functions)
colnames(bx.df) <- gsub("h(", "H(", colnames(bx.df), fixed=TRUE)
lm.bx    <- lm(survived ~ ., data=bx.df)
set.seed(2020)
earth.bx <- earth(survived ~ ., data=bx.df, degree=2)
printf("\nsummary(earth.bx):\n")
print(summary(earth.bx))
printf("\nevimp(earth.bx):\n")
print(evimp(earth.bx))
plot(earth.bx, info=TRUE)
plotmo(lm.bx) # Warning: Cannot determine which variables to plot in degree2 plots
plotmo(lm.bx, all2=TRUE, SHOWCALL=TRUE)
plotmo(earth.bx, pmethod="partdep", trace=2)

printf("\n=== put spaces into the column names of bx (for both response and predictors) ===\n")
spaced.bx <- bx.df
colnames(spaced.bx) <- gsub("-", " - ", colnames(spaced.bx), fixed=TRUE)
colnames(spaced.bx)[colnames(spaced.bx) == "survived"] <- "Survived = YES"
printf("\nhead(spaced.bx):\n")
print(head(spaced.bx))

lm.spaced.bx    <- lm(`Survived = YES` ~ ., data=spaced.bx)

set.seed(2020)
earth.spaced.bx <- earth(`Survived = YES` ~ ., data=spaced.bx, degree=2, trace=.5,
                         nfold=4, ncross=3, varmod.method="lm", pmethod="cv")
printf("\nsummary(earth.spaced.bx):\n")
print(summary(earth.spaced.bx))
printf("\nevimp(earth.spaced.bx):\n")
print(evimp(earth.spaced.bx))

set.seed(2020)
earth.glm.spaced.bx <- earth(`Survived = YES` ~ ., data=spaced.bx, degree=2, trace=.5,
                         glm=list(family="binomial"),
                         nfold=4, ncross=3, varmod.method="lm", pmethod="cv")
printf("\nsummary(earth.glm.spaced.bx):\n")
print(summary(earth.glm.spaced.bx))
printf("\nevimp(earth.glm.spaced.bx):\n")
print(evimp(earth.glm.spaced.bx))

options(warn=2)
expect.err(try(plotmo(lm.spaced.bx)),
           "Cannot determine which variables to plot in degree2 plots")
options(warn=1)

plotmo(lm.spaced.bx, do.par=2, SHOWCALL=TRUE)
plotres(lm.spaced.bx, do.par=0, which=c(1, 3))
par(org.par)

plotmo(earth.spaced.bx, degree1="sexmale", do.par=2, level=.8, SHOWCALL=TRUE)
plot(earth.spaced.bx, do.par=0, which=c(1, 3), info=TRUE, level=.8, type="earth")
par(org.par)

plot(earth.spaced.bx, versus="b:", info=TRUE, level=.8, type="earth", SHOWCALL=TRUE)

# following should be the same as previous page (since type="earth")
plotmo(earth.glm.spaced.bx, degree1="sexmale", do.par=2, level=.8, type="earth", SHOWCALL=TRUE)
plot(earth.glm.spaced.bx, do.par=0, which=1, info=TRUE, level=.8, type="earth")
# $$ TODO Following shouldn't cause Warning: Internal inconsistency: p$fit - fitted != 0
#         No warning if don't use glm=list(family="binomial") in call to earth
options(warn=2)
expect.err(try(plot(earth.glm.spaced.bx, do.par=0, which=3, info=TRUE, level=.8, type="earth")),
           "Internal inconsistency: p$fit != fitted")
options(warn=1)
plot(earth.glm.spaced.bx, do.par=0, which=3, info=TRUE, level=.8, type="earth")
par(org.par)

expect.err(try(plotmo(earth.glm.spaced.bx, level=.8)),
           "predict.earth: with earth-glm models, use type=\"earth\" when using the interval argument")

plotmo(earth.glm.spaced.bx, degree1="sexmale", do.par=2, SHOWCALL=TRUE)
plot(earth.glm.spaced.bx, do.par=0, which=c(1, 3), info=TRUE)
par(org.par)

printf("\n=== test combinations of variables in formula ===\n")

vdata <- data.frame(
    resp = 1:13,
    bool = c(F, F, F, F, F, T, T, T, T, T, T, T, T),
    ord  = ordered(c("ORD1", "ORD1", "ORD1",
                     "ORD1", "ORD1", "ORD1",
                     "ORD3", "ORD3", "ORD3",
                     "ORD2", "ORD2", "ORD2", "ORD2"),
                   levels=c("ORD1", "ORD3", "ORD2")),
    fac  = as.factor(c("FAC1", "FAC1", "FAC1",
                       "FAC2", "FAC2", "FAC2",
                       "FAC3", "FAC3", "FAC3",
                       "FAC1", "FAC2", "FAC3", "FAC3")),
    str  = c("STR1", "STR1", "STR1", # WILL BE TREATED LIKE A FACTOR
             "STR2", "STR2", "STR2",
             "STR3", "STR3", "STR3",
             "STR3", "STR3", "STR3", "STR3"),
    num  = c(1, 3, 2, 3, 4, 5, 6, 4, 5, 6.5, 3, 6, 5), # 7 unique values (but one is non integral)
    sqrt_num  = sqrt(c(1, 3, 2, 3, 4, 5, 6, 4, 5, 6.5, 3, 6, 5)),
    int  = c(1L, 1L, 3L, 3L, 4L, 4L, 3L, 5L, 3L, 6L, 7L, 8L, 10L), # 8 unique values
    date = as.Date(
           c("2018-08-01", "2018-08-02", "2018-08-03",
             "2018-08-04", "2018-08-05", "2018-08-06",
             "2018-08-07", "2018-08-08", "2018-08-08",
             "2018-08-08", "2018-08-10", "2018-08-11", "2018-08-11")),
    date_num = as.numeric(as.Date(
           c("2018-08-01", "2018-08-02", "2018-08-03",
             "2018-08-04", "2018-08-05", "2018-08-06",
             "2018-08-07", "2018-08-08", "2018-08-08",
             "2018-08-08", "2018-08-10", "2018-08-11", "2018-08-11"))))

vdata$off <- (1:nrow(vdata)) / nrow(vdata)

resp2 <- 13:1

vweights <- rep(1, length.out=nrow(vdata))
vweights[1] <- 2

set.seed(2020)
lognum.bool.ord.off <- earth(resp ~ log(num) + bool + ord + offset(off), degree=2, weights=vweights,
           data=vdata, pmethod="none", varmod.method="lm",
           nfold=2, ncross=3,
           trace=1)

printf("summary(lognum.bool.ord.off)\n")
print(summary(lognum.bool.ord.off))
cat("\nevimp(lognum.bool.ord.off)\n")
print(evimp(lognum.bool.ord.off))

plotmo(lognum.bool.ord.off, do.par=2, level=.8, SHOWCALL=TRUE)
plot(lognum.bool.ord.off, which=1, do.par=0)
par(org.par)

num.fac.sqrt.num.ord.bool <- earth(resp ~ num + int + fac + offset(off) + sqrt(num) + ord:bool - int,
           data=vdata, pmethod="none", trace=1)
plotmo(num.fac.sqrt.num.ord.bool, SHOWCALL=TRUE)
cat("\nevimp(num.fac.sqrt.num.ord.bool)\n")
print(evimp(num.fac.sqrt.num.ord.bool))

printf("\n=== unusual formulas, compare to lm ===\n")

lm1 <- lm(resp~ord+sqrt(as.numeric(fac)) + num+sqrt(num / 2)+I(2 * int)+date, data = vdata)

# same formula terms  as lm1 but in different order
earth1 <- earth(resp~sqrt(as.numeric(fac)) + ord + date + num + sqrt(.5 * num)+I(int / .5),
                data = vdata, linpreds=TRUE, thresh=0, penalty=-1)
cat("\nevimp(earth1)\n")
print(evimp(earth1))
plotmo(lm1,    SHOWCALL=TRUE)
plotmo(earth1, SHOWCALL=TRUE)
stopifnot(max(abs(sort(lm1$coef) - sort(earth1$coef))) < 1e-10)
stopifnot((summary(lm1)$r.squared - earth1$rsq) < 1e-10)
stopifnot(max(abs(predict(lm1, newdata=vdata[5,,drop=FALSE]) - predict(earth1, newdata=vdata[5,,drop=FALSE]))) < 1e-10)

fac.sqrt    <- earth(resp~sqrt(num)+fac, data = vdata, linpreds=TRUE, thresh=0, penalty=-1)
fac.sqrt_   <- earth(resp~sqrt_num+fac,  data = vdata, linpreds=TRUE, thresh=0, penalty=-1)
cat("\nevimp(fac.sqrt)\n")
print(evimp(fac.sqrt))
cat("\nevimp(fac.sqrt_)\n")
print(evimp(fac.sqrt_))
# as.vector to remove names (which are slightly different: sqrt(num) vs sqrt_num
stopifnot(identical(as.vector(fac.sqrt$coef), as.vector(fac.sqrt_$coef)))

newdata.extra <- vdata[3:5,] # extra variables unused in the model
newdata.extra$extra <- sqrt(newdata.extra[,1])
cat("\ncolnames(newdata.extra):", paste(colnames(newdata.extra)), "\n")

newd  <- vdata[3:5,c("num", "fac")]             # only variables used in the formula model
newd_ <- vdata[3:5,c("num", "sqrt_num", "fac")] # only variables used in the xy model

stopifnot(identical(predict(fac.sqrt, newdata=newdata.extra), predict(fac.sqrt_, newdata=newd_)))
stopifnot(identical(predict(fac.sqrt, newdata=newd), predict(fac.sqrt_, newdata=newd_)))
stopifnot(identical(predict(fac.sqrt, newdata=newd), predict(fac.sqrt_, newdata=newd_)))
stopifnot(identical(predict(fac.sqrt, newdata=newd), predict(fac.sqrt_, newdata=newd_)))

stopifnot(max(abs(predict(fac.sqrt, newdata=newdata.extra) - predict(fac.sqrt_, newdata=newdata.extra))) < 1e-10)
stopifnot(max(abs(predict(fac.sqrt, newdata=newdata.extra) - predict(fac.sqrt_, newdata=newdata.extra))) < 1e-10)

printf("\n=== two response model ===\n")

vdata.2resp <- vdata
resp2 <- 13:1
vdata.2resp$resp2 <- resp2

earth.2resp <- earth(resp+resp2~num+sqrt(num), data=vdata.2resp, weights=vweights, trace=1,
                     linpreds=TRUE, thresh=0, penalty=-1)
printf("\nsummary(earth.2resp)\n")
print(summary(earth.2resp))
cat("\nevimp(earth.2resp)\n")
print(evimp(earth.2resp))
par(mfrow = c(2, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0), oma=c(0,0,4,0))
# for formula models, plotmo plots a sinle plot for the effect of num
plotmo(earth.2resp, nresp=1, do.par=0, main="earth.2resp    nresp1")
title <- paste(
    "two-response model: resp+resp2~num+sqrt(num)\n",
    "the top row is for earth.formula models: the combined effect of num and sqrt(num) is plotted together\n",
    "the bottom row is for an earth.default model: num and sqrt(num) are plotted separately")
title(title, outer=TRUE, cex=.6)
plotmo(earth.2resp, nresp=2, do.par=0, main="earth.2resp    nresp2")

# put two response data mats into matrix form for earth.default and for lm
xmat <- vdata[,c("num", "sqrt_num"), drop=FALSE]
colnames(xmat) <- c("num", "sqrt(num)")
xmat <- as.matrix(xmat)
ymat <- vdata[, "resp", drop=FALSE]
ymat$resp2 <- resp2
ymat <- as.matrix(ymat)
earthxy.2resp <- earth(xmat, ymat, weights=vweights, trace=1,
                       linpreds=TRUE, thresh=0, penalty=-1)
printf("\nsummary(earthxy.2resp)\n")
print(summary(earthxy.2resp))
cat("\nevimp(earthxy.2resp)\n")
print(evimp(earthxy.2resp))
# for xy models, plotmo plots a separate plots for the effect of num and sqrt(num)
plotmo(earthxy.2resp, nresp=1, do.par=0)
# plotmo(earthxy.2resp, nresp=2, do.par=0)
stopifnot(identical(earth.2resp$coeff, earthxy.2resp$coeff))

lm.2resp <- lm(ymat~xmat, weights=vweights)
printf("\nsummary(lm.2resp)\n")
print(summary(lm.2resp))
options(warn=2) # treat warnings as errors
expect.err(try(plotmo(lm.2resp, nresp=1)),
           "the variable on the right side of the formula is a matrix or data.frame")
options(warn=1) # print warnings as they occur

# check that lm and earth coeffs are the same
# need order() below because coeffs appear in different row order in the coeff mat
earth.2resp.order <- order(earth.2resp$coeff[,1])
lm.order <- order(lm.2resp$coeff[,1])
stopifnot(max(abs(earth.2resp$coeff[earth.2resp.order] - lm.2resp$coeff[lm.order])) < 1e-10)

printf("\n=== test glm() with spaced.bx ===\n")

# glm requires response to be a factor (or two columns)
spaced.bx.fac <- spaced.bx
spaced.bx.fac$`surv fac` <- factor(ifelse(spaced.bx$`Survived = YES`, "yes", "no"), levels = c("yes", "no"))
spaced.bx.fac$`Survived = YES` <- NULL
glm.spaced.bx <- glm(`surv fac` ~ ., data=spaced.bx.fac, family="binomial")
printf("summary(glm.spaced.bx):\n")
print(summary(glm.spaced.bx))
plotmo(glm.spaced.bx, do.par=2)
plotres(glm.spaced.bx, which=3, do.par=0, info=TRUE, main="plotres(glm.spaced.bx,which=3")
# TODO why is Residuals-Vs-Fitted plot different for plotres and plot for glm models?
plot(glm.spaced.bx, which=1, caption="plot(glm.spaced.bx, which=1)")
par(org.par)
plotmo(glm.spaced.bx, all2=TRUE, degree2=c("sexmale", "pclass"), SHOW.CALL=TRUE, do.par=2)
plotmo(glm.spaced.bx, degree1=0, all2=TRUE, degree2=c("sexmale", "age"), do.par=0)
par(org.par) # TODO I think plot(glm.spaced.bx) doesn't restore the graphics params?

printf("\n=== test formulas which have a rhs variable which a matrix ===\n")
# This also tests that earth's naming of variables is the same as lm for such rhs variables
#
# TODO plotmo fails when rhs variable is a matrix --- would be nice to fix that

x_ <- etitanic[,"age",drop=FALSE]
x_$pclass <- etitanic$pclass
x_$pclass <- as.numeric(etitanic$pclass)
x_ <- as.matrix(x_)
y_ <- as.matrix(as.numeric(etitanic[,"survived"]))

earthxy.rhs.mat <- earth(x_, y_, degree=2, trace=1)
print(summary(earthxy.rhs.mat))
cat("\nevimp(earthxy.rhs.mat)\n")
print(evimp(earthxy.rhs.mat))
cat("\nearthxy.rhs.mat$modvars\n")
print(earthxy.rhs.mat$modvars)
plotmo(earthxy.rhs.mat, SHOWCALL=TRUE) # ok

earth.rhs.mat <- earth(y_ ~ x_, degree=2, trace=1)
print(summary(earth.rhs.mat))
cat("\nevimp(earth.rhs.mat)\n")
print(evimp(earth.rhs.mat))
cat("\nearth.rhs.mat$modvars\n")
print(earth.rhs.mat$modvars)
stopifnot(max(abs(earthxy.rhs.mat$coeff - earth.rhs.mat$coeff)) < 1e-15)
expect.err(try(plotmo(earth.rhs.mat)), # Warning: the variable on the right side of the formula is a matrix or data.frame
           "model.frame.default could not interpret the data passed to get.earth.x from model.matrix.earth from predict.earth")
expect.err(try(plotmo(earth.rhs.mat, all1=TRUE)), # still fails
           "model.frame.default could not interpret the data passed to get.earth.x from model.matrix.earth from predict.earth")

lm.rhs.mat <- lm(y_ ~ x_)
print(summary(lm.rhs.mat))
expect.err(try(plotmo(lm.rhs.mat)), # Warning: the variable on the right side of the formula is a matrix or data.frame
           "predict returned the wrong length (got 1046 but expected 50)")
expect.err(try(plotmo(lm.rhs.mat, all1=TRUE)), # still fails
           "predict returned the wrong length (got 1046 but expected 50)")
earth1.rhs.mat <- earth(y_ ~ x_, linpreds=TRUE, thresh=0, penalty=-1) # degree1
cat("\nevimp(earth1.rhs.mat)\n")
print(evimp(earth1.rhs.mat))
options(warn=2)
expect.err(try(plotmo(earth.rhs.mat)),
           "the variable on the right side of the formula is a matrix or data.frame")
expect.err(try(plotmo(earth.rhs.mat, all1=TRUE)), # still fails
           "the variable on the right side of the formula is a matrix or data.frame")
options(warn=1)
stopifnot(max(abs(sort(lm.rhs.mat$coeff) - sort(earth1.rhs.mat$coeff))) < 1e-12)
stopifnot(sort(rownames(lm.rhs.mat$coeff)) == sort(rownames(earth1.rhs.mat$coeff)))

x_nonames <- x_
colnames(x_nonames) <- NULL
lm.rhs.nonames <- lm(y_ ~ x_nonames)
print(summary(lm.rhs.nonames))
expect.err(try(plotmo(lm.rhs.nonames)), # Warning: the variable on the right side of the formula is a matrix or data.frame
           "predict returned the wrong length (got 1046 but expected 50)")
expect.err(try(plotmo(lm.rhs.nonames, all1=TRUE)), # still fails
           "predict returned the wrong length (got 1046 but expected 50)")
earth1.rhs.nonames <- earth(y_ ~ x_nonames, linpreds=TRUE, thresh=0, penalty=-1) # degree1
print(summary(earth1.rhs.nonames))
cat("\nevimp(earth1.rhs.nonames)\n")
print(evimp(earth1.rhs.nonames))
options(warn=2)
expect.err(try(plotmo(earth1.rhs.nonames)), # Warning: the variable on the right side of the formula is a matrix or data.frame
           "the variable on the right side of the formula is a matrix or data.frame")
expect.err(try(plotmo(earth1.rhs.nonames, all1=TRUE)), # still fails
           "the variable on the right side of the formula is a matrix or data.frame")
options(warn=1)
stopifnot(max(abs(sort(lm.rhs.nonames$coeff) - sort(earth1.rhs.nonames$coeff))) < 1e-12)
stopifnot(sort(rownames(lm.rhs.nonames$coeff)) == sort(rownames(earth1.rhs.nonames$coeff)))

printf("\n=== test handling consecutive '-' in formula ===\n")

options(warn=2)
lm.consec.minus <- lm(Volume~.--Girth, data=trees) # note double --
expect.err(try(plotmo(lm.consec.minus)),
           "Consecutive '-' in formula may cause problems")
earth.consec.minus <- earth(Volume~.--Girth, data=trees) # note double --
cat("\nsummary(earth.consec.minus)\n")
print(summary(earth.consec.minus))
cat("\nevimp(earth.consec.minus)\n")
print(evimp(earth.consec.minus))
expect.err(try(plotmo(earth.consec.minus)),
           "Consecutive '-' in formula may cause problems")
options(warn=1)

printf("\n=== test rpart() with spaced.bx ===\n")

library(rpart.plot)
rpart.mod <- rpart(`Survived = YES` ~ ., data=spaced.bx)
printf("\nprint(rpart.rules(rpart.mod))\n")
print(rpart.rules(rpart.mod))
set.seed(2020)
plotmo(rpart.mod, do.par=2, degree1=c("sexmale", "pclass3rd"), degree2=2, pt.col="red")
plotres(rpart.mod, do.par=0, which=c(1,3))
par(org.par)

printf("\n=== tibble, class \"Date\", and ndiscrete ===\n")

library(tibble)
library(lubridate)
tib1 <- tibble(y    = c(1, 1, 2, 3), # even number of variables
               bool = c(F, F, F, T),
               date = c(ymd('2018-08-01'), ymd('2018-08-02'), ymd('2018-08-03'),
                        ymd('2018-08-03')))
cat("class tib1$date: ", class(tib1$date), "\n")
mod.tib1 <- lm(y ~ ., data = tib1)
plotmo(mod.tib1, col.response=2, all2=TRUE, ticktype="d", do.par=2, caption="mod.tib1: Dates ndiscrete=default 5")
plotmo(mod.tib1, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45)
par(org.par)

plotmo(mod.tib1, col.response=2, all2=TRUE, ticktype="d", do.par=2, ndiscrete=2, caption="mod.tib1: Dates ndiscrete=2")
plotmo(mod.tib1, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45, ndiscrete=2)
par(org.par)

plotmo(mod.tib1, col.response=2, all2=TRUE, ticktype="d", do.par=2, ndiscr=1, caption="mod.tib1: Dates ndiscrete=1")
plotmo(mod.tib1, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45, ndiscrete=2)
par(org.par)

tib2 <- tibble(y    = c(1, 1, 2, 3, 4), # odd number of variables
               bool = c(F, F, F, T, T),
               date = c(ymd('2018-08-01'), ymd('2018-08-02'), ymd('2018-08-03'),
                        ymd('2018-08-03'), ymd('2018-08-04')))
mod.tib2 <- lm(y ~ ., data = tib2)
plotmo(mod.tib2, col.response=2, all2=TRUE, ticktype="d", do.par=2, caption="mod.tib2: Dates ndiscrete=default 5")
plotmo(mod.tib2, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45)
par(org.par)

plotmo(mod.tib2, col.response=2, all2=TRUE, ticktype="d", do.par=2, ndiscrete=2, caption="mod.tib2: Dates ndiscrete=2")
plotmo(mod.tib2, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45, ndiscrete=2)
par(org.par)

plotmo(mod.tib2, col.response=2, all2=TRUE, ticktype="d", do.par=2, ndiscr=1, caption="mod.tib2: Dates ndiscrete=1")
plotmo(mod.tib2, col.response=2, degree1=0, all2=TRUE, ticktype="d", do.par=0, theta=-45, ndiscrete=2)
par(org.par)

source("test.epilog.R")
