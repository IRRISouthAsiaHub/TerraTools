## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ----echo=FALSE,eval=TRUE,warning=FALSE, message=FALSE------------------------
library(spdep)

## ----echo=TRUE,eval=TRUE------------------------------------------------------
library(spdep)
nc <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
#st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)

## ----echo=TRUE----------------------------------------------------------------
nc$ft.SID74 <- sqrt(1000)*(sqrt(nc$SID74/nc$BIR74) + sqrt((nc$SID74+1)/nc$BIR74))
nc$both <- factor(paste(nc$L_id, nc$M_id, sep=":"))

## ----echo=FALSE---------------------------------------------------------------
is_tmap <- FALSE
if (require(tmap, quietly=TRUE)) is_tmap <- TRUE
tmap4 <- FALSE
if(is_tmap) tmap4 <- packageVersion("tmap") >= "3.99"

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
gal_file <- system.file("weights/ncCC89.gal", package="spData")[1]
ncCC89 <- read.gal(gal_file, region.id=nc$FIPSNO)

## ----echo=TRUE----------------------------------------------------------------
sids.nhbr30.dist <- nbdists(ncCC89, cbind(nc$east, nc$north))
sids.nhbr <- listw2sn(nb2listw(ncCC89, glist=sids.nhbr30.dist, style="B", zero.policy=TRUE))
dij <- sids.nhbr[,3]
n <- nc$BIR74
el1 <- min(dij)/dij
el2 <- sqrt(n[sids.nhbr$to]/n[sids.nhbr$from])
sids.nhbr$weights <- el1*el2
if (packageVersion("spdep") >= "1.3.1") {
  sids.nhbr.listw <- sn2listw(sids.nhbr, style="B", zero.policy=TRUE)
} else {
  sids.nhbr.listw <- sn2listw(sids.nhbr)
}

## ----echo=TRUE----------------------------------------------------------------
nc$ft.NWBIR74 <- sqrt(1000)*(sqrt(nc$NWBIR74/nc$BIR74) + sqrt((nc$NWBIR74+1)/nc$BIR74))

## ----echo=TRUE----------------------------------------------------------------
lm_nc <- lm(ft.SID74 ~ 1, data=nc)
outl <- which.max(rstandard(lm_nc))
as.character(nc$NAME[outl])

## ----echo=TRUE----------------------------------------------------------------
W <- listw2mat(sids.nhbr.listw)
W.4 <- W[-outl, -outl]
sids.nhbr.listw.4 <- mat2listw(W.4)
nc2 <- nc[!(1:length(nc$CNTY_ID) %in% outl),]

## ----echo=TRUE----------------------------------------------------------------
library(spatialreg)
ecarIaw <- spautolm(ft.SID74 ~ 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIaw)

## ----echo=TRUE----------------------------------------------------------------
ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIIaw)

## ----echo=TRUE----------------------------------------------------------------
ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIVaw)

## ----eval=is_tmap, echo=TRUE--------------------------------------------------
nc2$fitIV <- fitted.values(ecarIVaw)
if (tmap4) {
  tm_shape(nc2) + tm_polygons(fill="fitIV", fill.scale=tm_scale(values="brewer.yl_or_br"), fill.legend=tm_legend(position=tm_pos_in("left", "bottom"), frame=FALSE, item.r = 0), lwd=0.01)
} else {
tm_shape(nc2) + tm_fill("fitIV")
}

## ----echo=TRUE----------------------------------------------------------------
ecarIawll <- spautolm(ft.SID74 ~ 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR", llprof=seq(-0.1, 0.9020532358, length.out=100))
plot(ll ~ lambda, ecarIawll$llprof, type="l")

