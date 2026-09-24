# create README.html from README.md
# the paths below assume that this file is in the plotmo/inst/slowtests directory
library(rpart.plot)
library(rmarkdown)
rmarkdown::render("../../README.md", output_dir="../../.#")
if(!interactive())
    q(runLast=FALSE)
