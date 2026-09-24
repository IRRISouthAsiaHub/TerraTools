### R code from vignette source 'GPA1guide.Rnw'

###################################################
### code chunk number 1: GPA1guide.Rnw:22-24
###################################################
options(continue="  ")
pdf.options(pointsize = 8)


###################################################
### code chunk number 2: GPA1guide.Rnw:45-46
###################################################
library("GPArotation")


###################################################
### code chunk number 3: GPA1guide.Rnw:90-101
###################################################
data("Harman", package = "GPArotation")

# Calling a rotation directly
qHarman <- quartimax(Harman8)

# Equivalently, via the wrapper function
qHarman <- GPFRSorth(Harman8, method = "quartimax")

# Two equivalent ways to access the rotated loadings
loadings(qHarman)       # via extractor function (recommended)
qHarman$loadings        # via direct list access


###################################################
### code chunk number 4: GPA1guide.Rnw:127-138
###################################################
data("CCAI", package = "GPArotation")
fa.un <- factanal(factors = 3, covmat = CCAI_R, n.obs = 461,
                  rotation = "none")

res.quart <- quartimax(fa.un)
max(abs(res.quart$loadings %*% t(res.quart$Th) - fa.un$loadings))

res.obli <- oblimin(fa.un, normalize = TRUE, randomStarts = 15)
max(abs(res.obli$loadings %*% t(res.obli$Th) - fa.un$loadings))

max(abs(res.obli$loadings - fa.un$loadings %*% solve(t(res.obli$Th))))


###################################################
### code chunk number 5: GPA1guide.Rnw:146-147
###################################################
max(abs(res.obli$Phi - t(res.obli$Th) %*% res.obli$Th))


###################################################
### code chunk number 6: GPA1guide.Rnw:166-177
###################################################
set.seed(334)
res <- quartimin(Harman8, normalize = TRUE, randomStarts = 100)

# Raw unsorted loadings
loadings(res)

# Sorted loadings via print (factors reordered and signs adjusted)
res.sorted <- print(res)

# Once sorted, repeated calls to print are stable
max(abs(print(res.sorted)$loadings - res.sorted$loadings)) == 0  # TRUE


###################################################
### code chunk number 7: GPA1guide.Rnw:216-221
###################################################
res.obli <- oblimin(Harman8, normalize = TRUE, randomStarts = 100)
# Pattern matrix (unsorted)
summary(res.obli, Structure = FALSE)
# Structure matrix
summary(res.obli, Structure = TRUE)


###################################################
### code chunk number 8: GPA1guide.Rnw:264-284
###################################################
plotPatternStructure <- function(pattern, structure,
                                  labels = NULL,
                                  main = "Pattern vs Structure") {
  k   <- ncol(pattern)
  col <- palette.colors(k, palette = "Okabe-Ito")
  par(mfrow = c(1, k))
  for (j in 1:k) {
    lims <- range(c(pattern[, j], structure[, j]))
    plot(pattern[, j], structure[, j],
         xlim = lims, ylim = lims,
         xlab = "Pattern loading",
         ylab = "Structure loading",
         main = paste(if (!is.null(labels)) labels[j]
                      else paste("Factor", j)),
         pch  = 19, col = col[j])
    abline(0, 1, lty = 2, col = "grey60")
    abline(h = 0, col = "grey80")
    abline(v = 0, col = "grey80")
  }
}


###################################################
### code chunk number 9: GPA1guide.Rnw:288-294
###################################################
res.obli.s  <- print(res.obli)
Pattern     <- loadings(res.obli.s)
Structure   <- Pattern %*% res.obli.s$Phi

plotPatternStructure(Pattern, Structure,
                     labels = c("Factor 1", "Factor 2"))


###################################################
### code chunk number 10: GPA1guide.Rnw:323-367
###################################################
plotRotationLandscape <- function(A, method = "quartimax", n = 20000,
                                   main = NULL, ...) {
  if (ncol(A) != 2)
    stop("plotRotationLandscape only works for 2-factor solutions.")

  vgQfun_fn <- get(paste("vgQ", method, sep = "."),
                   envir = asNamespace("GPArotation"))

  if (is.null(main))
    main <- paste("Rotation landscape:", method)

  theta  <- seq(0, 2 * pi, length.out = n)
  f_vals <- numeric(n)

  for (i in seq_along(theta)) {
    Tmat      <- matrix(c( cos(theta[i]), sin(theta[i]),
                          -sin(theta[i]), cos(theta[i])), 2, 2)
    L         <- A %*% Tmat
    VgQ       <- do.call(vgQfun_fn, append(list(L), list(...)))
    f_vals[i] <- VgQ$f
  }

  min_idx <- which.min(f_vals)

  plot(theta, f_vals,
       type  = "l",
       lwd   = 2,
       main  = main,
       xlab  = expression(theta ~ "(radians)"),
       ylab  = "f",
       xaxt  = "n")
  axis(1, at     = c(0, pi/2, pi, 3*pi/2, 2*pi),
           labels = c("0", expression(pi/2), expression(pi),
                      expression(3*pi/2), expression(2*pi)))
  abline(h   = f_vals[min_idx], col = "grey80", lty = 2)
  points(theta[min_idx], f_vals[min_idx],
         col = "tomato", pch = 19, cex = 1.5)
  text(theta[min_idx], f_vals[min_idx],
       labels = paste0("min f = ", round(f_vals[min_idx], 4)),
       pos    = 4,
       cex    = 0.8)

  invisible(data.frame(theta = theta, f = f_vals))
}


###################################################
### code chunk number 11: GPA1guide.Rnw:373-374
###################################################
data(Harman, package = "GPArotation")


###################################################
### code chunk number 12: GPA1guide.Rnw:378-383
###################################################
par(mfrow = c(2, 2))
plotRotationLandscape(Harman8, method = "quartimax")
plotRotationLandscape(Harman8, method = "varimax")
plotRotationLandscape(Harman8, method = "bentler")
plotRotationLandscape(Harman8, method = "entropy")


###################################################
### code chunk number 13: GPA1guide.Rnw:412-434
###################################################
res <- plotRotationLandscape(Harman8, method = "simplimax", k = 4)

# Find all local minima by sign changes in the discrete derivative
df         <- diff(res$f)
local_mins <- which(df[-length(df)] < 0 & df[-1] > 0)

# Classify minima: red = symmetry-equivalent global, blue = genuine local
f_rounded  <- round(res$f[local_mins], 4)
f_min      <- min(f_rounded)
is_global  <- f_rounded == f_min
dot_cols   <- ifelse(is_global, "tomato", "steelblue")
n_global   <- sum(is_global)
n_local    <- sum(!is_global)

points(res$theta[local_mins], res$f[local_mins],
       col = dot_cols, pch = 19, cex = 0.9)

legend("topright",
       legend = c(paste0("global minimum (x", n_global, ")"),
                  paste0("local minima (x",   n_local,  ")")),
       col    = c("tomato", "steelblue"),
       pch    = 19, bty = "n", cex = 0.8)


###################################################
### code chunk number 14: GPA1guide.Rnw:438-441
###################################################
cat("Total minima found:              ", length(local_mins), "\n")
cat("Symmetry-equivalent global min:  ", n_global, "\n")
cat("Genuine local minima:            ", n_local,  "\n")


###################################################
### code chunk number 15: GPA1guide.Rnw:473-504
###################################################
data(Harman, package = "GPArotation")
res.quart   <- quartimax(Harman8)
res.oblimin <- oblimin(Harman8)

L.quart   <- abs(loadings(res.quart))
L.oblimin <- abs(loadings(res.oblimin))

ord.quart   <- order(L.quart[, 1],   decreasing = TRUE)
ord.oblimin <- order(L.oblimin[, 1], decreasing = TRUE)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

barplot(t(L.quart[ord.quart, ]),
        beside      = TRUE,
        ylim        = c(0, 1),
        main        = "Quartimax",
        ylab        = "Absolute loading",
        xlab        = "Variable (sorted by Factor 1)",
        legend.text = c("Factor 1", "Factor 2"),
        args.legend = list(x = "topright"),
        col         = c("steelblue", "tomato"))

barplot(t(L.oblimin[ord.oblimin, ]),
        beside      = TRUE,
        ylim        = c(0, 1),
        main        = "Oblimin",
        ylab        = "Absolute loading",
        xlab        = "Variable (sorted by Factor 1)",
        legend.text = c("Factor 1", "Factor 2"),
        args.legend = list(x = "topright"),
        col         = c("steelblue", "tomato"))


###################################################
### code chunk number 16: GPA1guide.Rnw:525-577
###################################################
plotSortedLoadings <- function(..., labels = NULL, col = NULL,
                                main = "Sorted Absolute Loadings",
                                ylab = "Absolute loading",
                                xlab = "Rank") {
  solutions <- list(...)

  for (i in seq_along(solutions)) {
    if (!inherits(solutions[[i]], "GPArotation"))
      stop("Argument ", i, " is not a GPArotation object.")
  }

  n <- length(solutions)

  if (is.null(labels))
    labels <- paste("Solution", seq_len(n))
  if (is.null(col))
    col <- palette.colors(n, palette = "Okabe-Ito")

  sorted_loadings <- lapply(solutions, function(x)
    sort(abs(as.vector(x$loadings)), decreasing = FALSE))

  all_values <- unlist(sorted_loadings)
  max_len    <- max(sapply(sorted_loadings, length))

  plot(NULL,
       xlim = c(1, max_len),
       ylim = c(0, max(all_values)),
       main = main,
       xlab = xlab,
       ylab = ylab,
       las  = 1)

  abline(h = seq(0, 1, by = 0.1), col = "grey90", lty = 1)

  for (i in seq_len(n)) {
    lines(seq_along(sorted_loadings[[i]]), sorted_loadings[[i]],
          col = col[i], lwd = 2)
    points(seq_along(sorted_loadings[[i]]), sorted_loadings[[i]],
           col = col[i], pch = 19, cex = 0.6)
  }

  legend("topleft", legend = labels, col = col, lwd = 2, pch = 19,
         bty = "n")

  invisible(sorted_loadings)
}

# Example
data(Harman, package = "GPArotation")
res.quart   <- quartimax(Harman8)
res.oblimin <- oblimin(Harman8)
res.geomin  <- geominT(Harman8)


###################################################
### code chunk number 17: GPA1guide.Rnw:583-585
###################################################
plotSortedLoadings(res.quart, res.oblimin, res.geomin,
                   labels = c("Quartimax", "Oblimin", "Geomin"))


###################################################
### code chunk number 18: GPA1guide.Rnw:642-673
###################################################
origdigits <- options("digits")
options(digits = 2)
trBritain <- matrix(c(.783,-.163,.811,.202,.724,.209,.850,.064,
  -.031,.592,-.028,.723,.388,.434,.141,.808,.215,.709),
  byrow = TRUE, ncol = 2)
trGermany <- matrix(c(.778,-.066,.875,.081,.751,.079,.739,.092,
  .195,.574,-.030,.807,-.135,.717,.125,.738,.060,.691),
  byrow = TRUE, ncol = 2)
# orthogonal rotation of trGermany towards trBritain
trx <- targetT(trGermany, Target = trBritain)
# Factor loadings after target rotation
trx
# Differences between loadings matrices after rotation
y <- trx$loadings - trBritain
print(y, digits = 1)
# Square root of the mean squared difference per item
sqrt(apply((y^2), 1, mean))
# Square root of the mean squared difference per factor
sqrt(apply((y^2), 2, mean))
# Identity coefficient per factor after rotation
2 * colSums(trx$loadings * trBritain) /
  (colSums(trx$loadings^2) + colSums(trBritain^2))
# Additivity coefficient per factor after rotation
diag(2 * cov(trx$loadings, trBritain)) /
  diag(var(trx$loadings) + var(trBritain))
# Proportionality coefficient per factor after rotation
colSums(trBritain * trx$loadings) /
  sqrt(colSums(trBritain^2) * colSums(trx$loadings^2))
# Correlation for each factor after rotation
diag(cor(trBritain, trx$loadings))
options(digits = origdigits$digits)


###################################################
### code chunk number 19: GPA1guide.Rnw:685-708
###################################################
plot(trBritain[, 1], trBritain[, 2],
     xlim = c(-0.3, 1.0), ylim = c(-0.3, 1.0),
     xlab = "Factor 1", ylab = "Factor 2",
     main = "Target Rotation: Germany towards Britain",
     pch = 19, col = "steelblue", cex = 1.2)
abline(h = 0, lty = 2, col = "grey70")
abline(v = 0, lty = 2, col = "grey70")
points(trGermany[, 1], trGermany[, 2],
       pch = 17, col = "tomato", cex = 1.2)
points(loadings(trx)[, 1], loadings(trx)[, 2],
       pch = 15, col = "orange", cex = 1.2)
for (i in 1:nrow(trGermany)) {
  arrows(trGermany[i, 1], trGermany[i, 2],
         loadings(trx)[i, 1], loadings(trx)[i, 2],
         length = 0.08, col = "grey60")
}
legend("topright",
       legend = c("Britain (varimax rotated)",
                  "East Germany (varimax rotated)",
                  "East Germany (rotated towards Britain)"),
       col    = c("steelblue", "tomato", "orange"),
       pch    = c(19, 17, 15),
       bty    = "n")


###################################################
### code chunk number 20: GPA1guide.Rnw:735-749
###################################################
A <- matrix(c(.664, .688, .492, .837, .705, .82, .661, .457, .765, .322,
  .248, .304, -0.291, -0.314, -0.377, .397, .294, .428, -0.075, .192, .224,
  .037, .155, -.104, .077, -.488, .009), ncol = 3)
# using targetT
SPA <- matrix(c(rep(NA, 6), .7, .0, .7, rep(0, 3), rep(NA, 7),
  0, 0, NA, 0, rep(NA, 4)), ncol = 3)
xt <- targetT(A, Target = SPA)
# using pstT
SPApst <- matrix(c(rep(0, 6), .7, .0, .7, rep(0, 3), rep(0, 7),
  0, 0, 0, 0, rep(0, 4)), ncol = 3)
SPAW <- matrix(c(rep(0, 6), rep(1, 6), rep(0, 7), 1, 1, 0, 1,
  rep(0, 4)), ncol = 3)
xpst <- pstT(A, Target = SPApst, W = SPAW)
max(abs(loadings(xt) - loadings(xpst)))


###################################################
### code chunk number 21: GPA1guide.Rnw:766-770
###################################################
data("CCAI", package = "GPArotation")
factanal(factors = 3, covmat = CCAI_R, n.obs = 461, rotation = "infomaxT")
factanal(factors = 3, covmat = CCAI_R, n.obs = 461, rotation = "infomaxT",
  control = list(rotate = list(normalize = TRUE, eps = 1e-6)))


###################################################
### code chunk number 22: GPA1guide.Rnw:782-786
###################################################
data("CCAI", package = "GPArotation")
fa.un <- factanal(factors = 3, covmat = CCAI_R, n.obs = 461,
                  rotation = "none")
oblimin(fa.un, randomStarts = 100)


###################################################
### code chunk number 23: GPA1guide.Rnw:799-818
###################################################
data("WansbeekMeijer", package = "GPArotation")
fa.un <- factanal(factors = 3, covmat = NetherlandsTV,
                  normalize = TRUE, rotation = "none")

set.seed(42)
res.cf <- cfQ(fa.un, kappa = 0.3, normalize = TRUE, randomStarts = 100)
res.cf

if (getRversion() >= "4.5.1") {
  set.seed(42)
  fa.cf <- factanal(factors = 3, covmat = NetherlandsTV,
                    normalize = TRUE, rotation = "cfQ",
                    control = list(rotate = list(kappa = 0.3,
                                                 randomStarts = 100)))
  cat("Maximum difference in loadings:\n")
  print(max(abs(abs(res.cf$loadings) - abs(fa.cf$loadings))))
} else {
  cat("Use the two-step procedure above for correct results.\n")
}


###################################################
### code chunk number 24: browne-table5
###################################################
Lambda <- matrix(c(
    0.90, 0.00, 0.00,
    0.90, 0.00, 0.00,
    0.90, 0.00, 0.00,
    0.90, 0.00, 0.00,
    0.00, 0.80, 0.00,
    0.00, 0.80, 0.00,
    0.00, 0.80, 0.00,
    0.00, 0.00, 0.70,
    0.00, 0.00, 0.70,
    0.00, 0.00, 0.70,
    0.10, 0.10, 0.00,
    0.10, 0.00, 0.10
), nrow = 12, ncol = 3, byrow = TRUE,
   dimnames = list(paste0("Var", 1:12), paste0("Fac", 1:3)))

h          <- sqrt(rowSums(Lambda^2))
kaiser_wts <- round(1 / h, 2)
al_cm      <- GPArotation:::NormalizingWeight(Lambda,
                                              normalize = "CM")[, 1]
cm_wts     <- round(pmin(pmax(h / al_cm, 0), 1), 2)

print(cbind(round(Lambda, 2),
            "K Wts."  = kaiser_wts,
            "CM Wts." = cm_wts))


###################################################
### code chunk number 25: browne-cm-varimax
###################################################
rot_cm <- GPForth(Lambda, method = "varimax", normalize = "CM")
print(round(unclass(rot_cm$loadings), 2))


###################################################
### code chunk number 26: browne-table6
###################################################
rhos     <- c(0.0, 0.3, 0.5)
patterns <- list()

for (rho in rhos) {
    Phi_true        <- matrix(rho, 3, 3); diag(Phi_true) <- 1
    R_common        <- Lambda %*% Phi_true %*% t(Lambda)
    ev              <- eigen(R_common)
    A               <- ev$vectors[, 1:3] %*% diag(sqrt(ev$values[1:3]))
    rownames(A)     <- paste0("Var", 1:12)
    rot             <- GPFoblq(A, method = "cf",
                               methodArgs = list(kappa = 1/12),
                               normalize  = "CM")
    P               <- unclass(rot$loadings)
    col_order       <- order(colSums(P^2), decreasing = TRUE)
    P               <- P[, col_order]
    for (j in 1:3) if (sum(P[, j]) < 0) P[, j] <- -P[, j]
    colnames(P)     <- paste0("Fac", 1:3)
    patterns[[as.character(rho)]] <- round(P, 2)
}

t6 <- cbind(patterns[["0"]], patterns[["0.3"]], patterns[["0.5"]])
colnames(t6) <- c("rho=0: F1", "F2", "F3",
                  "rho=.3: F1", "F2", "F3",
                  "rho=.5: F1", "F2", "F3")
print(t6)
