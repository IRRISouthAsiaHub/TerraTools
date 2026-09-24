### R code from vignette source 'GPA4fitstats.Rnw'

###################################################
### code chunk number 1: GPA4fitstats.Rnw:49-52
###################################################
library("GPArotation")
data("GriffithMulaik", package = "GPArotation")
n.obs <- 523


###################################################
### code chunk number 2: GPA4fitstats.Rnw:64-73
###################################################
for (k in 3:8) {
  fa.k  <- factanal(factors = k, covmat = GriffithMulaik,
                    n.obs = n.obs, rotation = "none")
  res.k <- quartimax(fa.k)
  fit.k <- GPArotation:::calc_fitstats(res.k)
  cat(sprintf("k = %d  RMSEA = %.4f  [%.4f, %.4f]  SRMR = %.4f\n",
              k, fit.k$RMSEA, fit.k$RMSEA.l, fit.k$RMSEA.u,
              fit.k$SRMR))
}


###################################################
### code chunk number 3: GPA4fitstats.Rnw:139-156
###################################################

fa.un <- factanal(factors = 6, covmat = GriffithMulaik,
                  n.obs = n.obs, rotation = "none")
res.vm <- Varimax(fa.un)
res.ob <- oblimin(fa.un, randomStarts = 100)
res.t1 <- tandemI(fa.un)

fit.vm <- GPArotation:::calc_fitstats(res.vm)
fit.ob <- GPArotation:::calc_fitstats(res.ob)
fit.t1 <- GPArotation:::calc_fitstats(res.t1)

cat(sprintf("Varimax  SRMR = %.4f  RMSEA = %.4f\n",
            fit.vm$SRMR, fit.vm$RMSEA))
cat(sprintf("Oblimin  SRMR = %.4f  RMSEA = %.4f\n",
            fit.ob$SRMR, fit.ob$RMSEA))
cat(sprintf("TandemI  SRMR = %.4f  RMSEA = %.4f\n",
            fit.t1$SRMR, fit.t1$RMSEA))


###################################################
### code chunk number 4: GPA4fitstats.Rnw:166-167
###################################################
GPArotation:::audit_residuals(res.ob)


###################################################
### code chunk number 5: GPA4fitstats.Rnw:178-179
###################################################
plot(res.ob, "residuals")


###################################################
### code chunk number 6: GPA4fitstats.Rnw:205-206
###################################################
print(res.ob)


###################################################
### code chunk number 7: GPA4fitstats.Rnw:215-216
###################################################
summary(res.ob)


###################################################
### code chunk number 8: GPA4fitstats.Rnw:229-245
###################################################
res.gm <- geominQ(fa.un, randomStarts = 100)

auc.ob <- GPArotation:::calc_AUC(res.ob)$AUC_mean
auc.vm <- GPArotation:::calc_AUC(res.vm)$AUC_mean
auc.gm <- GPArotation:::calc_AUC(res.gm)$AUC_mean

hp.ob  <- GPArotation:::calc_hyperplane(res.ob)$HP_pct
hp.vm  <- GPArotation:::calc_hyperplane(res.vm)$HP_pct
hp.gm  <- GPArotation:::calc_hyperplane(res.gm)$HP_pct

cat(sprintf("%-12s  AUC = %.3f  Hyperplane = %.1f pct\n",
            "Varimax",  auc.vm, hp.vm))
cat(sprintf("%-12s  AUC = %.3f  Hyperplane = %.1f pct\n",
            "Oblimin",  auc.ob, hp.ob))
cat(sprintf("%-12s  AUC = %.3f  Hyperplane = %.1f pct\n",
            "GeominQ",  auc.gm, hp.gm))


###################################################
### code chunk number 9: GPA4fitstats.Rnw:256-263
###################################################
sim.un <- GPArotation:::calc_simplicity(fa.un)
sim.ob <- GPArotation:::calc_simplicity(res.ob)

cat(sprintf("%-12s  Hoffman = %.3f  Gini = %.3f  Bentler = %.3f\n",
            "Unrotated", sim.un$Hoffman, sim.un$Gini, sim.un$Bentler))
cat(sprintf("%-12s  Hoffman = %.3f  Gini = %.3f  Bentler = %.3f\n",
            "Oblimin",   sim.ob$Hoffman, sim.ob$Gini, sim.ob$Bentler))
