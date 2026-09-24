# test.parsnip.R: test the parsnip package with earth and other models
# Stephen Milborrow Sep 2020 Petaluma

source("test.prolog.R")
options(warn=1) # print warnings as they occur
library(earth)
cat("loading parsnip libraries\n") # these libraries take several seconds to load
library(tidymodels, quietly=TRUE, verbose=FALSE)
library(timetk)
library(lubridate)
cat("loaded parsnip libraries\n")
cat("parsnip version:", as.character(packageVersion("parsnip")[[1]]), "\n")

vdata <- data.frame(
    resp    = 1:23,
    bool = c(F, F, F, F, F, T, T, T, T, T, T, T, T, F, F, T, T, T, T, T, T, T, T),
    ord  = ordered(c("ORD1", "ORD1", "ORD1",
                     "ORD1", "ORD1", "ORD1",
                     "ORD1", "ORD3", "ORD1",
                     "ORD2", "ORD2", "ORD2", "ORD2",
                     "ORD2", "ORD2", "ORD2",
                     "ORD3", "ORD3", "ORD3",
                     "ORD2", "ORD2", "ORD2", "ORD2"),
                   levels=c("ORD1", "ORD3", "ORD2")),
    fac  = as.factor(c("FAC1", "FAC1", "FAC1",
                       "FAC2", "FAC2", "FAC2",
                       "FAC3", "FAC1", "FAC1",
                       "FAC1", "FAC2", "FAC2", "FAC2",
                       "FAC2", "FAC2", "FAC2",
                       "FAC3", "FAC3", "FAC3",
                       "FAC1", "FAC3", "FAC3", "FAC3")),
    str  = c("STR1", "STR1", "STR1", # WILL BE TREATED LIKE A FACTOR
             "STR1", "STR1", "STR1",
             "STR2", "STR2", "STR2",
             "STR3", "STR3", "STR2", "STR3",
             "STR2", "STR3", "STR2",
             "STR3", "STR3", "STR3",
             "STR3", "STR3", "STR3", "STR3"),
    num  = c(1, 9, 2, 3, 14, 5, 6, 4, 5, 6.5, 3, 6, 5,
             3, 4, 5, 6, 4, 5, 16.5, 3, 16, 15),
    sqrt_num  = sqrt(
           c(1, 9, 2, 3, 14, 5, 6, 4, 5, 6.5, 3, 6, 5,
             3, 4, 5, 6, 4, 5, 16.5, 3, 16, 15)),
    int  = c(1L, 1L, 3L, 3L, 4L, 4L, 3L, 5L, 3L, 6L, 7L, 8L, 10L,
             13L, 14L, 3L, 13L, 5L, 13L, 16L, 17L, 18L, 11L),
    date = as.Date(
           c("2018-08-01", "2018-08-02", "2018-08-03",
             "2018-08-04", "2018-08-05", "2018-08-06",
             "2018-08-07", "2018-08-08", "2018-08-08",
             "2018-08-10", "2018-08-10", "2018-08-11", "2018-08-11",
             "2018-08-11", "2018-08-12", "2018-08-13",
             "2018-08-10", "2018-08-15", "2018-08-17",
             "2018-08-04", "2018-08-19", "2018-08-03", "2018-08-18")),
    date_num = as.numeric(as.Date(
           c("2018-08-01", "2018-08-02", "2018-08-03",
             "2018-08-04", "2018-08-05", "2018-08-06",
             "2018-08-07", "2018-08-08", "2018-08-08",
             "2018-08-10", "2018-08-10", "2018-08-11", "2018-08-11",
             "2018-08-11", "2018-08-12", "2018-08-13",
             "2018-08-10", "2018-08-15", "2018-08-17",
             "2018-08-04", "2018-08-19", "2018-08-03", "2018-08-18"))))

set.seed(2020)
splits <- initial_time_split(vdata, prop=.9)

#--- lm ----------------------------------------------------------------------

lm1 <- lm(resp~num+fac:int+date+ord+str, data=training(splits))
cat("lm1:\n")
print(summary(lm1))
set.seed(2020)
lmpar <- linear_reg(mode = "regression") %>%
         set_engine("lm") %>%
         fit(resp~num+fac:int+date+ord+str, data = training(splits))
stopifnot(identical(lm1$coeff, lmpar$fit$coeff))

predict.lm1   <- predict(lm1, testing(splits))
predict.lmpar <- lmpar %>% predict(testing(splits))
stopifnot(all(predict.lm1 == predict.lmpar))

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
plotmo(lm1,        do.par=2, SHOWCALL=TRUE)
plotres(lm1,       which=c(3,1), do.par=FALSE)
plotmo(lmpar,      do.par=2, SHOWCALL=TRUE)
plotres(lmpar,     which=c(3,1), do.par=FALSE)
plotmo(lmpar$fit,  do.par=2, SHOWCALL=TRUE)
plotres(lmpar$fit, which=c(3,1), do.par=FALSE)
par(org.par)

lmpar.sqrtnum <- linear_reg(mode = "regression") %>%
         set_engine("lm") %>%
         fit(resp~sqrt(num), data = training(splits))
#$$ TODO
# expect.err(try(plotmo(lmpar.sqrtnum)),
#            "cannot get the original model predictors")

#--- earth -------------------------------------------------------------------

# note that sqrt(num) is ok, unlike parsnip models for lm and rpart
earth1 <- earth(resp~sqrt(num)+int+ord:bool+fac+str+date, degree=2,
                data=training(splits), pmethod="none")
cat("earth1:\n")
print(summary(earth1))
set.seed(2020)
earthpar <- mars(mode = "regression", prune_method="none", prod_degree=2) %>%
            set_engine("earth") %>%
            fit(resp~sqrt(num)+int+ord:bool+fac+str+date, data = training(splits))
cat("earthpar:\n")
print(earthpar)
cat("summary(earthpar$fit)\n")
print(summary(earthpar$fit))
stopifnot(identical(earth1$coeff, earthpar$fit$coeff))

predict.earth1 <- predict(earth1, testing(splits))
predict.earthpar <- earthpar %>% predict(testing(splits))
stopifnot(all(predict.earth1 == predict.earthpar))

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
plotmo(earth1, do.par=2, pt.col=3, SHOWCALL=TRUE)
set.seed(2020)
plotres(earth1, which=c(1,3), do.par=FALSE, pt.col=3, legend.pos="topleft")
par(org.par)

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
plotmo(earthpar, do.par=2, pt.col=3, SHOWCALL=TRUE)
set.seed(2020)
plotres(earthpar, which=c(1,3), do.par=FALSE, pt.col=3, legend.pos="topleft")
par(org.par)

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
plotmo(earthpar$fit, do.par=2, pt.col=3, SHOWCALL=TRUE)
set.seed(2020)
plotres(earthpar$fit, which=c(1,3), do.par=FALSE, pt.col=3, legend.pos="topleft")
par(org.par)

#--- rpart -------------------------------------------------------------------

library(rpart)
library(rpart.plot)
rpart1 <- rpart(resp~num+fac+int+date+ord+str, data=training(splits),
                control=rpart.control(minsplit=1, cp=.0001))
cat("\nrpart.rules(rpart1)\n")
print(rpart.rules(rpart1))

set.seed(2020)
# TODO note need of model=TRUE below (needed only for further processing with e.g. plotmo)
rpartpar <- decision_tree(mode = "regression", min_n=1, cost_complexity=.0001) %>%
             set_engine("rpart", model=TRUE) %>%
             fit(resp~num+fac+int+date+ord+str, data = training(splits))
cat("\nrpart.rules(rpartpar$fit)\n")
print(rpart.rules(rpartpar$fit))

predict.rpart1   <- predict(rpart1, testing(splits))
predict.rpartpar <- rpartpar %>% predict(testing(splits))
stopifnot(all(predict.rpart1 == predict.rpartpar))

par(mfrow = c(3, 3), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
plotmo(rpart1,        do.par=2, SHOWCALL=TRUE, trace=0)
plotres(rpart1,       which=c(3,1), do.par=FALSE)
plotmo(rpartpar,      do.par=2, SHOWCALL=TRUE, trace=0)
plotres(rpartpar,     which=c(3,1), do.par=FALSE)
plotmo(rpartpar$fit,  do.par=2, SHOWCALL=TRUE)
plotres(rpartpar$fit, which=c(3,1), do.par=FALSE)
par(org.par)

# TODO note that this differs from the above rpart model in that we don't use model=TRUE
rpartpar.nosavemodel <- decision_tree(mode = "regression", min_n=1, cost_complexity=.0001) %>%
             set_engine("rpart") %>%
             fit(resp~num+fac+int+date+str, data = training(splits))

cat("\nrpart.rules(rpartpar.nosavemodel$fit)\n")
options(warn=2)
expect.err(try(rpart.rules(rpartpar.nosavemodel$fit)),
           "Cannot retrieve the data used to build the model")
options(warn=1)
expect.err(try(plotmo(rpartpar.nosavemodel)),
           "Cannot plot parsnip rpart model: need model=TRUE in call to rpart")

rpart.sqrtnum <- decision_tree(mode = "regression", min_n=1, cost_complexity=.0001) %>%
             set_engine("rpart", model=TRUE) %>%
             fit(resp~sqrt(num)+fac+int+date+ord+str, data = training(splits))
cat("\nrpart.rules(rpart.sqrtnum$fit)\n")
print(rpart.rules(rpart.sqrtnum$fit)) # ok
#$$ TODO
# expect.err(try(plotmo(rpart.sqrtnum)),
#            "cannot get the original model predictors")

#-----------------------------------------------------------------------------------
# Test fix for github bug report https://github.com/tidymodels/parsnip/issues/341
# (fixed Sep 2020)

cat("===m750a first example===\n")
set.seed(2020)
m750a <- m4_monthly %>%
    filter(id == "M750") %>%
    select(-id)
print(m750a) # a tibble
set.seed(2020)
splits_a <- initial_time_split(m750a, prop = 0.9)
earth_m750a <- earth(log(value) ~ as.numeric(date) + month(date, label = TRUE), data = training(splits_a), degree=2)
print(summary(earth_m750a))
set.seed(2020)
model_m750a <- mars(mode = "regression", prod_degree=2) %>%
    set_engine("earth") %>%
    fit(log(value) ~ as.numeric(date) + month(date, label = TRUE), data = training(splits_a))
print(summary(model_m750a$fit))
stopifnot(identical(earth_m750a$coeff, model_m750a$fit$coeff))
predict_earth_m750a <- predict(earth_m750a, newdata=testing(splits_a)[1:3,])
predict_m750a <- model_m750a %>% predict(testing(splits_a)[1:3,])
stopifnot(max(c(9.238049628, 9.240535151, 9.232361834) - predict_m750a) < 1e-8)
stopifnot(max(predict_earth_m750a - predict_m750a) < 1e-20)

par(mfrow = c(2, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
set.seed(2020)
plotmo(model_m750a, trace=2, do.par=FALSE, pt.col="green", main="model_m750a", SHOWCALL=TRUE)
set.seed(2020)
plotmo(model_m750a$fit, trace=1, do.par=FALSE, pt.col="green", main="model_m750a$fit", SHOWCALL=TRUE)
set.seed(2020)
plotmo(earth_m750a, trace=1, do.par=FALSE, pt.col="green", main="earth_m750a", SHOWCALL=TRUE)
par(org.par)

cat("===m750a second example===\n")
set.seed(2020)
m750b <- m4_monthly %>%
    filter(id == "M750") %>%
    select(-id) %>%
    rename(date2 = date)
print(m750b) # tibble
set.seed(2020)
splits_b <- initial_time_split(m750b, prop = 0.9)
set.seed(2020)
model_m750b <- mars(mode = "regression") %>%
    set_engine("earth") %>%
    fit(log(value) ~ as.numeric(date2) + month(date2, label = TRUE), data = training(splits_b))
# new data that only contains the feature "date" as a predictor
future_data <- m750b %>% future_frame(date2, .length_out = "3 years")
print(future_data) # a tibble with a single column of class "Date"
stopifnot(class(future_data[,1,drop=TRUE]) == "Date")
predict_m750a <- model_m750b %>% predict(new_data = future_data)

par(mfrow = c(2, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
set.seed(2020)
plotmo(model_m750b, trace=2, do.par=FALSE, pt.col="green", main="model_m750b", SHOWCALL=TRUE)
set.seed(2020)
plotmo(model_m750b$fit, trace=1, do.par=FALSE, pt.col="green", main="model_m750b$fit", SHOWCALL=TRUE)
par(org.par)

#-----------------------------------------------------------------------------------
# multiple response earth model

data(etitanic)

etit <- etitanic
etit$survived <- factor(ifelse(etitanic$survived == 1, "yes", "no"),
                        levels = c("yes", "no"))
etit$notsurvived <- factor(ifelse(etitanic$survived == 0, "notsurvived", "survived"),
                        levels = c("notsurvived", "survived"))
set.seed(2020)
earth_tworesp <- earth(survived + notsurvived ~ ., data=etit, degree=2)
print(summary(earth_tworesp))

# TODO following commented out because parsnip (version 0.1.5) says "'+' not meaningful for factors"
# set.seed(2020)
# mars_tworesp <- mars(mode = "regression", prod_degree=2) %>%
#          set_engine("earth") %>%
#          fit(survived + notsurvived~., data=etit)
# print(summary(mars_tworesp))
# print(summary(mars_tworesp$fit))
#
# stopifnot(identical(earth_tworesp$coeff, mars_tworesp$fit$coeff))
#
# predict.earth_tworesp <- predict(earth_tworesp, etit[3:6,])
# predict.mars_tworesp <-  mars_tworesp %>% predict(etit[3:6,])
# stopifnot(all(predict.earth_tworesp == predict.mars_tworesp))
#
# plotmo(earth_tworesp, trace=0, nresponse=1, SHOWCALL=TRUE)
# plotmo(mars_tworesp, trace=0, nresponse=1, SHOWCALL=TRUE)
# plotmo(mars_tworesp, trace=0, nresponse=2, SHOWCALL=TRUE)

source("test.epilog.R")
