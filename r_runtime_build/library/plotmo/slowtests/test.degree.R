# test.pre.R: test the degree1 and degree2 and related args

source("test.prolog.R")
library(earth)
library(plotmo)

# test character degree1 and degree2 (added in plotmo version 1.3-0)

data(ozone1)
a80 <- earth(O3~., data=ozone1, degree=2)
plotmo(a80, degree1="i", degree2="t",
       caption='degree1="i", degree2="t"')
plotmo(a80, degree1="^temp$", degree2="^dpg$",
       caption='degree1="^temp$", degree2="^dpg$"')
# Expect Warning: "nonesuch1" in degree1 does not regex-match any variables, ditto for degree2
plotmo(a80, degree1=c("temp", "nonesuch1"), degree2="vis",
       caption='degree1=c("temp", "nonesuch1"), degree2="vis")')
# Expect above warnings and also Warning: nothing to plot
plotmo(a80, degree1="nonesuch1", degree2="nonesuch2")

# tests for plotmo version 3.3.7 (degree1 and degree2 handling changed)

data(etitanic)
a81 <- earth(survived~., data=etitanic, degree=2)
options(warn=1) # print warnings as they occur
plotmo(a81)

# degree1 tests
par(mfrow=c(3,3), mar=c(1,2.5,2,1), oma=c(0,0,4,0))
plotmo(a81, do.par=FALSE, degree1="pclass", degree2=0, main='degree1="pclass"',
       caption="test degree1 with strings")
options(warn=2) # treat warnings as errors
expect.err(try(plotmo(a81, do.par=FALSE, degree1="survived", degree2=0)), '"survived" in degree1 does not regex-match any names')
options(warn=1) # print warnings as they occur
plotmo(a81, do.par=FALSE, degree1="sibsp", degree2=0, main='degree1="sibsp"')
# parch does not appear in the standard degree1 plotmo plots, but we can still specify it explictly
plotmo(a81, do.par=FALSE, degree1="parch", degree2=0, trace=0, main='degree1="parch"')
plotmo(a81, do.par=FALSE, degree1=c("sibsp", "pclass"), degree2=0, main='degree1=c("sibsp", "pclass")')
par(org.par)

# degree2 tests
par(mfrow=c(3,3), mar=c(1,2.5,2,1), oma=c(0,0,4,0))
plotmo(a81, do.par=FALSE, degree1=0, degree2="pclass", main='degree2="pclass"',
       caption="test degree2 with two strings")
plotmo(a81, do.par=FALSE, degree1=0, degree2=c("age", "se"), persp.theta=-35,
       main='degree2=c("age", "se")\npersp.theta=-35')
plotmo(a81, do.par=FALSE, degree1=0, degree2="ag", main='degree2="ag"')
plotmo(a81, do.par=FALSE, degree1=0, degree2=c("sex", "sibsp"), main='degree2=c("sex", "sibsp"')
plotmo(a81, do.par=FALSE, degree1=0, degree2=c("sibsp", "sex"), main='degree2=c("sibsp", "sex")')
options(warn=2) # treat warnings as errors
expect.err(try(plotmo(a81, do.par=FALSE, degree1=0, degree2=c("pclass", "nonesuch"))), "\"nonesuch\" in degree2 does not regex-match any names")
expect.err(try(plotmo(a81, do.par=FALSE, degree1=0, degree2=c("nonesuch1", "nonesuch2"))), "\"nonesuch1\" in degree2 does not regex-match any names")
expect.err(try(plotmo(a81, do.par=FALSE, degree1=0, degree2=c("nonesuch", "pclass"))), "\"nonesuch\" in degree2 does not regex-match any names")
options(warn=1) # print warnings as they occur
par(org.par)

par(mfrow=c(2,2), mar=c(1,2.5,2,1), oma=c(0,0,4,0))

# check that order of strings in two string degree2 is observed
cat('\n\ndegree2=c("age", "se"):\n')
plotmo(a81, do.par=FALSE, degree1=0,
       degree2=c("age", "se"), main='degree2=c("age", "se")')
cat('\n\ndegree2=c("se", "age"):\n')
plotmo(a81, do.par=FALSE, degree1=0,
       degree2=c("se", "age"), main='degree2=c("se", "age")')

# check handling of bad strings in two string degree2
cat('\n\ndegree2=c("nonesuch", "age"):\n')
try(plotmo(a81, do.par=FALSE, degree1=0,
    degree2=c("nonesuch", "age"), main='degree2=c("nonesuch", "age")'))
cat('\n\ndegree2=c("age", "nonesuch"):\n')
try(plotmo(a81, do.par=FALSE, degree1=0,
       degree2=c("age", "nonesuch"),
       main='degree2=c("age", "nonesuch")'))
cat('\n\ndegree2=c("nevermore", "nonesuch"):\n')
try(plotmo(a81, do.par=FALSE, degree1=0,
       degree2=c("nevermore", "nonesuch"),
       main='degree2=c("nevermore", "nonesuch")'))
# follow should still plot the degree1 plot even though degree2 spec is wrong
cat('\n\ndegree1=1, degree2=c("nevermore", "nonesuch"):\n')
try(plotmo(a81, do.par=FALSE, degree1=1,
       degree2=c("nevermore", "nonesuch"),
       main='degree1=1\ndegree2=c("nevermore", "nonesuch")'))

# expect warning: both elements of degree2 are the same
cat('\n\ndegree2=c("sex", "sex"):\n')
try(plotmo(a81, do.par=FALSE, degree1=0,
       degree2=c("sex", "sex"),
       main='degree1=1\ndegree2=c("sex", "sex")'))

par(org.par)

source("test.epilog.R")
