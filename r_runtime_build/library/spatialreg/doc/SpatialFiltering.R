## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ----echo=TRUE----------------------------------------------------------------
library(spdep)
require("sf", quietly=TRUE)
if (packageVersion("spData") >= "2.3.2") {
    NY8 <- sf::st_read(system.file("shapes/NY8_utm18.gpkg", package="spData"))
} else {
    NY8 <- sf::st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData"))
    sf::st_crs(NY8) <- "EPSG:32618"
    NY8$Cases <- NY8$TRACTCAS
}
NY_nb <- read.gal(system.file("weights/NY_nb.gal", package="spData"), override.id=TRUE)

## ----echo=TRUE----------------------------------------------------------------
library(spatialreg)
nySFE <- SpatialFiltering(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, nb=NY_nb, style="W", verbose=FALSE)
nylmSFE <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+fitted(nySFE), data=NY8)
summary(nylmSFE)
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
anova(nylm, nylmSFE)

## ----echo=TRUE----------------------------------------------------------------
NYlistwW <- nb2listw(NY_nb, style = "W")
set.seed(111)

## ----echo=TRUE----------------------------------------------------------------
nyME <- ME(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, offset=log(POP8), family="poisson", listw=NYlistwW, alpha=0.46)

## ----echo=TRUE----------------------------------------------------------------
nyME
NY8$eigen_1 <- fitted(nyME)[,1]
NY8$eigen_2 <- fitted(nyME)[,2]

## ----echo=TRUE----------------------------------------------------------------
#gry <- brewer.pal(9, "Greys")[-1]
plot(NY8[,c("eigen_1", "eigen_2")])

## ----echo=TRUE----------------------------------------------------------------
nyglmME <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+fitted(nyME), data=NY8, family="poisson")
summary(nyglmME)
nyGLMp <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8,family="poisson")
anova(nyGLMp, nyglmME, test="Chisq")

