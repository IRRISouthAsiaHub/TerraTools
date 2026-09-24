# limod.methods.R: Additional method functions for the linmod example.
#
# See www.milbo.org/doc/modguide.pdf.
# This software may be freely used.

variable.names.linmod <- function(object, ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    names(coef(object))
}
case.names.linmod <- function(object, ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    names(residuals(object))
}
nobs.linmod <- function(object, use.fall.back = FALSE, ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    NROW(object$residuals)
}
deviance.linmod <- function(object, ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    sum(residuals(object)^2)
}
model.frame.linmod <- function(formula, ...)
{
    stopifnot(inherits(formula, "linmod"))
    if(is.null(formula$terms)) # model built with linmod.default?
        stop("model.frame cannot be used on linmod models built without a formula")
    else
        model.frame.default(formula, ...)
}
model.matrix.linmod <- function(object, data = NULL, ...)
{
    stopifnot(inherits(object, "linmod"))
    if(is.null(data))
        data <- model.frame.linmod(object)
    model.matrix.default(object, data = data, ...)
}
logLik.linmod <- function(object, REML = FALSE, ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    stopifnot(!REML) # linmod does not save qr hence cannot do REML
    res <- object$residuals
    p <- object$rank
    n <- length(res)
    w <- rep.int(1, n)
    n0 <- n
    val <- .5* (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) +
                                   log(sum(w*res^2))))
    attr(val, "nall") <- n0
    attr(val, "nobs") <- n
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    val
}
estfun.linmod <- function (x, ...) # for sandwich package
{
    stopifnot(inherits(x, "linmod"))
    stop.if.dot.arg.used(...)
    xmat <- model.matrix(x)
    res <- residuals(x)
    rval <- as.vector(res) * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}
plot.linmod <- function(x, main = NULL, ...) # dots are passed to plot()
{
    stopifnot(inherits(x, "linmod"))
    call.as.char <- paste0(deparse(x$call, control = NULL, nlines = 5),
                           sep = " ", collapse = " ")
    plot(fitted(x), residuals(x), xlab = "Fitted values", ylab = "Residuals",
         main = if(is.null(main)) substr(call.as.char, 1, 50) else main,
         ...)
    smooth <- lowess(fitted(x), residuals(x), f = .5)
    lines(smooth$x, smooth$y, col = 2)
}
