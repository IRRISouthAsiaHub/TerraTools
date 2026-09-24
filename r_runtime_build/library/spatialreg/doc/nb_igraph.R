## -----------------------------------------------------------------------------
library(spatialreg)

## -----------------------------------------------------------------------------
dothis <- TRUE
if (!suppressPackageStartupMessages(require(sf, quietly=TRUE))) {
  message("install the sf package")
  dothis <- FALSE
}
if (dothis) {
  sf_extSoftVersion()
}

## ----echo=dothis, eval=dothis-------------------------------------------------
library(sf)
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1])
row.names(columbus)[1:10]

## ----echo=dothis, eval=dothis-------------------------------------------------
nb_q <- spdep::poly2nb(columbus)
nb_q
attr(nb_q, "region.id")[1:10]
spdep::is.symmetric.nb(nb_q)

## ----echo=dothis, eval=dothis-------------------------------------------------
col2 <- spdep::droplinks(nb_q, 21)
nb_q[[21]]
col2[[21]]
col2
spdep::is.symmetric.nb(col2)
coords <- st_coordinates(st_centroid(st_geometry(columbus)))
plot(nb_q, coords, col="grey")
plot(col2, coords, add=TRUE)

## ----echo=dothis, eval=dothis-------------------------------------------------
nb_B <- spdep::nb2listw(col2, style="B", zero.policy=TRUE)
nb_B$style

## ----echo=dothis, eval=dothis-------------------------------------------------
library(spatialreg)
library(Matrix)
B <- as(nb_B, "CsparseMatrix")
all(B == t(B))
str(B)
rownames(B)[1:10]

## ----echo=dothis, eval=dothis-------------------------------------------------
nb_B1 <- spdep::mat2listw(as(as(B, "TsparseMatrix"), "CsparseMatrix"),
    style="B", zero.policy=TRUE)
nb_B1$style
all.equal(nb_B1$neighbours, col2, check.attributes=FALSE)
all.equal(attr(nb_B1$neighbours, "region.id"), attr(nb_B$neighbours, "region.id"))

## ----echo=dothis, eval=dothis-------------------------------------------------
rho <- 0.1
do_spatialreg <- FALSE
if (requireNamespace("spatialreg", quietly=TRUE)) do_spatialreg <- TRUE
if (do_spatialreg) sum(log(1 - rho * spatialreg::eigenw(nb_B)))

## ----echo=dothis, eval=dothis-------------------------------------------------
n <- nrow(B)
I <- Diagonal(n)
class(I - rho * B)
c(determinant(I - rho * B, logarithm=TRUE)$modulus)

## ----echo=dothis, eval=dothis-------------------------------------------------
nW <- -B
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho), sqrt=TRUE)$modulus))

## ----echo=dothis, eval=dothis-------------------------------------------------
nb_W <- spdep::nb2listw(col2, style="W", zero.policy=TRUE)
W <- as(nb_W, "CsparseMatrix")
str(W)
all(W == t(W))

## ----echo=dothis, eval=dothis-------------------------------------------------
set.seed(1)
x <- runif(n)
r1 <- as.numeric(W %*% x)
r2 <- lag(nb_W, x, zero.policy=TRUE)
all.equal(r1, r2, check.attributes=FALSE)
plot(x, r1, ylim=c(0,1))
c(x[21], r1[21])

## ----echo=dothis, eval=dothis-------------------------------------------------
rho <- 0.5
sum(log(1 - rho * eigenw(nb_W)))
class(I - rho * W)
c(determinant(I - rho * W, logarithm=TRUE)$modulus)

## ----echo=dothis, eval=dothis-------------------------------------------------
LU <- lu(I - rho * W)
sum(log(abs(diag(slot(LU, "U")))))

## ----echo=dothis, eval=dothis-------------------------------------------------
d <- attr(nb_W$weights, "comp")$d
all.equal(d, spdep::card(col2))

## ----echo=dothis, eval=dothis-------------------------------------------------
dW <- Diagonal(n, d) %*% W
all(dW == t(dW))
isd <- Diagonal(n, 1/sqrt(d))
isd[21,21]
Ws <- as(isd %*% dW %*% isd, "symmetricMatrix")
rowSums(Ws)[21]
class(Ws)
c(determinant(I - rho * Ws, logarithm=TRUE)$modulus)

## ----echo=dothis, eval=dothis-------------------------------------------------
1/range(spatialreg::eigenw(nb_B))
if (!require("RSpectra", quietly=TRUE)) dothis <- FALSE

## ----echo=dothis, eval=dothis-------------------------------------------------
1/c(eigs(B, k=1, which="SR")$values, eigs(B, k=1, which="LR")$values)

## ----echo=dothis, eval=dothis-------------------------------------------------
1/range(eigenw(nb_W))
1/Re(c(eigs(W, k=1, which="SR")$values, eigs(W, k=1, which="LR")$values))

## ----echo=dothis, eval=dothis-------------------------------------------------
class(B)
object.size(B)
if (!require("igraph", quietly=FALSE)) dothis <- FALSE

## ----echo=dothis, eval=dothis-------------------------------------------------
g1 <- graph_from_adjacency_matrix(B, mode="undirected")
class(g1)
object.size(g1)

## ----echo=dothis, eval=dothis-------------------------------------------------
# Matrix 1.4-2 vulnerability work-around
ow <- options("warn")$warn
options("warn"=2L)
B1 <- try(as_adjacency_matrix(g1), silent=TRUE)
if (!inherits(B1, "try-error")) print(class(B1))
if (!inherits(B1, "try-error")) print(object.size(B1))
if (!inherits(B1, "try-error")) print(all.equal(B, as(B1, "CsparseMatrix")))
options("warn"=ow)

## ----echo=dothis, eval=dothis-------------------------------------------------
res <- spdep::n.comp.nb(col2)
table(res$comp.id)

## ----echo=dothis, eval=dothis-------------------------------------------------
c1 <- components(g1)
c1$no == res$nc
all.equal(c1$membership, res$comp.id)
all.equal(c1$csize, c(table(res$comp.id)), check.attributes=FALSE)

## ----echo=dothis, eval=dothis-------------------------------------------------
W <- as(spdep::nb2listw(col2, style="W", zero.policy=TRUE), "CsparseMatrix")
g1W <- graph_from_adjacency_matrix(W, mode="directed", weighted="W")
c1W <- components(g1W)
all.equal(c1W$membership, res$comp.id)

## ----echo=dothis, eval=dothis-------------------------------------------------
is_connected(g1)
dg1 <- diameter(g1)
dg1
sp_mat <- distances(g1)
str(sp_mat)

## ----echo=dothis, eval=dothis-------------------------------------------------
nbl10 <- spdep::nblag(col2, maxlag=10)
vals <- sapply(nbl10, function(x) sum(spdep::card(x)))
zero <- which(vals == 0)
zero[1]-1

## ----echo=dothis, eval=dothis-------------------------------------------------
lmat <- lapply(nbl10[1:(zero[1]-1)], spdep::nb2mat, style="B", zero.policy=TRUE)
mat <- matrix(0, n, n)
for (i in seq(along=lmat)) mat = mat + i*lmat[[i]]
mat[mat==0] <- Inf
diag(mat) <- 0
all.equal(mat, sp_mat, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-------------------------------------------------
nb_r <- spdep::cell2nb(7, 7, type="rook")
nb_rW <- spdep::nb2listw(nb_r, style="W")
spdep:::find_q1_q2(nb_rW)

## ----echo=dothis, eval=dothis-------------------------------------------------
1/range(Re(eigenw(similar.listw(nb_rW))))

## ----echo=dothis, eval=dothis-------------------------------------------------
spdep:::find_q1_q2(nb_W)
1/range(Re(eigenw(similar.listw(nb_W))))

