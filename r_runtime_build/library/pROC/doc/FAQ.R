## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----legacy-axis-comparison, fig.width=10, fig.height=5, fig.align='center', out.width='100%', message=FALSE, warning=FALSE----
library(pROC)
data(aSAH)

# Create a ROC curve
rocobj <- roc(aSAH$outcome, aSAH$s100b)

# Plot with default axis (Specificity from 1 to 0)
par(mfrow = c(1, 2))
plot(rocobj, main = "Default: Specificity (1 to 0)")

# Plot with legacy axis (1-Specificity from 0 to 1)
plot(rocobj, main = "Legacy: 1-Specificity (0 to 1)", 
     legacy.axes = TRUE)
par(mfrow = c(1, 1))

## ----glm-example, fig.width=5, fig.height=5, fig.align='center', message=FALSE, warning=FALSE----
# Fit a logistic regression model
model <- glm(outcome ~ s100b + wfns, data = aSAH, family = binomial)

# Get predicted probabilities
predictions <- predict(model, type = "response")

# Create ROC curve with the predictions
roc_glm <- roc(aSAH$outcome, predictions)

# Plot and compute AUC
plot(roc_glm, main = "ROC curve from GLM predictions")
auc(roc_glm)

## ----default-empty-space, fig.width=7, fig.height=5, fig.align='center', message=FALSE, warning=FALSE----
rocobj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
plot(rocobj, main = "Default ROC plot (non-squared, with empty space)")

## ----remove-empty-space, fig.width=7, fig.height=5, fig.align='center', message=FALSE, warning=FALSE----
par(pty = "s")
rocobj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
plot(rocobj, main = "ROC curve with square plotting region")
par(pty = "m")  # Reset to default

## ----eval = FALSE-------------------------------------------------------------
# # For file output
# png("roc.png", width = 480, height = 480)
# par(mar = c(4, 4, 4, 4) + 0.1)
# plot(roc(aSAH$outcome, aSAH$s100b, percent = TRUE))
# dev.off()

## ----free-aspect-ratio, fig.width=7, fig.height=5, fig.align='center', message=FALSE, warning=FALSE----
rocobj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
plot(rocobj, main = "ROC curve with free aspect ratio", asp = NA)

## ----eval = FALSE-------------------------------------------------------------
# ci.coords(aSAH$outcome, aSAH$s100b, x="best", input = "threshold",
#           ret="threshold")
# # 95% CI (2000 stratified bootstrap replicates):
# #            2.5%   50%  97.5%
# # threshold 0.115 0.205 0.4854

## ----eval = FALSE-------------------------------------------------------------
# roc(outcome ~ ndka, data = aSAH)
# # Not enough distinct predictions to compute area under the ROC curve.
# #
# # Error in roc(outcome ~ ndka, data = aSAH) : unused argument (data = aSAH)

