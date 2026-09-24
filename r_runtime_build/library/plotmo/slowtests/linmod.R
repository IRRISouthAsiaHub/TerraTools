# linmod.R: Example S3 linear model.
#
# See www.milbo.org/doc/modguide.pdf.
# This software may be freely used.

linmod <- function(...) UseMethod("linmod")

linmod.default <- function(x = stop("no 'x' argument"),
                           y = stop("no 'y' argument"),
                           keep = FALSE,
                           ...)
{
    stop.if.dot.arg.used(...)
    xmat <- as.matrix(x)
    # use name "(Intercept)" here so coef names match linmod.formula
    x <- cbind("(Intercept)" = 1, xmat)
    fit <- linmod.fit(x, y)
    fit$call <- match.call()
    if(keep) {
        fit$x <- xmat
        # save y as a one-column matrix, so can use colname to save response name
        colname <- deparse(substitute(y))[1]
        colname <- gsub(" ", "", substr(colname, 1, 100)) # strip spaces, truncate
        fit$y <- as.matrix(y, ncol = 1)
        colnames(fit$y) <- colname
    }
    fit
}
linmod.formula <- function(formula = stop("no 'formula' argument"),
                           data = parent.frame(),
                           keep = FALSE,
                           ...)
{
    stop.if.dot.arg.used(...)
    if(is.matrix(data))             # allow data to be a matrix
        data <- as.data.frame(data) # will create colnames V1 V2 V3 if necessary
    # note that na.action=na.pass because we will catch NAs later
    # in linmod.fit, for uniformity with linmod.default
    mf <- model.frame(formula = formula, data = data, na.action = na.pass)
    terms <- attr(mf, "terms")
    x <- model.matrix(terms, mf)
    y <- model.response(mf)
    fit <- linmod.fit(x, y)
    fit$call <- match.call()
    fit$terms <- terms
    fit$xlevels <- .getXlevels(terms, mf) # for use by predict.linmod
    if(keep)
        fit$data <- data
    fit
}
linmod.fit <- function(x = stop("no 'x' argument"),
                       y = stop("no 'y' argument"),
                       ...)
{
    # internal function, not for the casual user
    # if model has an intercept, the first col of x must be intercept (all 1s)

    stop.if.dot.arg.used(...)
    x <- check.linmod.x(x)
    y <- check.linmod.y(x, y)
    fit <- do.linmod.fit(x, y)
    class(fit) <- "linmod"
    fit
}
check.linmod.x <- function(x)
{
    if(!is.matrix(x))
        stop("'x' is not a matrix or could not be converted to a matrix")
    if(NROW(x) == 0 || NCOL(x) == 0)
        stop("'x' is empty")
    if(anyNA(x))
        stop("NA in 'x'")
    # checking just the first column of x suffices because all columns
    # of a matrix have the same type
    # we allow is.logical because qr etc. treat logical vars as numeric
    if(!is.numeric(x[,1]) && !is.logical(x[,1]))
        stop("non-numeric column in 'x'")
    # ensure all columns in x are named (needed for names in vcov etc.)
    # use the same naming convention as lm (prefix for unnamed cols is "V")
    missing.colnames <-
        if(is.null(colnames(x))) 1:NCOL(x)
        else                     nchar(colnames(x)) == 0
    colnames(x)[missing.colnames] <-
        c("(Intercept)",
          paste("V", seq_len(NCOL(x) - 1), sep = ""))[missing.colnames]
    duplicated <- which(duplicated(colnames(x)))
    if(length(duplicated))
        stop("column name \"", colnames(x)[duplicated[1]],
             "\" in 'x' is duplicated")
    x
}
check.linmod.y <- function(x, y)
{
    # as.vector(as.matrix(y)) is necessary when y is a data.frame
    # (because as.vector alone on a data.frame returns a data.frame)
    y <- as.vector(as.matrix(y))
    if(length(y) == 0)
        stop("'y' is empty")
    if(anyNA(y))
        stop("NA in 'y'")
    if(!is.numeric(y) && !is.logical(y))
        stop("'y' is not numeric or logical")
    if(length(y) != nrow(x))
        stop("nrow(x) is ", nrow(x), " but length(y) is ", length(y))
    y
}
do.linmod.fit <- function(x, y)
{
    # workhorse function for fitting linear models
    # essential processing and sanity checks on x and y are already completed
    # x is a numeric matrix, y is a numeric vector

    qx <- qr(x)                         # QR-decomposition of x
    if(qx$rank < ncol(x))
        stop("'x' is singular (it has ", ncol(x),
             " columns but its rank is ", qx$rank, ")\n  colnames(x): ",
             paste0(colnames(x), collapse=' '))
    coef <- solve.qr(qx, y)             # compute (x'x)^(-1) x'y
    stopifnot(!anyNA(coef))             # NA impossible after rank check above
    df.residual <- max(0, nrow(x) - ncol(x)) # degrees of freedom
    sigma2 <- sum((y - x %*% coef)^2) / df.residual # variance of residuals
    vcov <- sigma2 * chol2inv(qx$qr)    # covar mat is sigma^2 * (x'x)^(-1)
    fitted.values <- qr.fitted(qx, y)

    colnames(vcov) <- rownames(vcov) <- colnames(x)
    names(fitted.values) <- rownames(x)
    colnames(coef) <- colnames(y)

    # returned fields match lm's fields
    list(coefficients  = coef,
         residuals     = y - fitted.values,
         rank          = qx$rank,
         fitted.values = fitted.values,
         vcov          = vcov,
         sigma         = sqrt(sigma2),
         df.residual   = df.residual)
}
predict.linmod <- function(object = stop("no 'object' argument"),
                           newdata = NULL,
                           type = "response",
                           ...)
{
    stopifnot(inherits(object, "linmod"))
    stop.if.dot.arg.used(...)
    match.arg(type, "response") # the type argument is not yet supported
    if(is.null(newdata))
        yhat <- fitted(object)
    else {
        if(NROW(newdata) == 0)
            stop("'newdata' is empty")      # preempt obscure message later

        x <- if(is.null(object$terms))      # model built with linmod.default?
                process.newdata(object, newdata)
            else                            # model built with linmod.formula
                process.newdata.formula(object, newdata)

        # The following tests suffice to catch all illegal input.  However
        # they aren't ideal in that they don't always direct you to the root
        # cause of the problem (i.e. the error messages aren't always optimal).
        nvar <- length(object$coefficients) - 1 # nbr vars, -1 for intercept
        if(ncol(x) - 1 != nvar)
            stop("ncol(newdata) is ", ncol(x) - 1, " but should be ", nvar)
        if(anyNA(x))
            stop("NA in 'newdata'")
        if(!is.numeric(x[,1]) && !is.logical(x[,1]))
            stop("non-numeric column in 'newdata' (after processing)")

        yhat <- as.vector(do.predict.linmod(object, x))
        names(yhat) <- rownames(x)
    }
    yhat
}
process.newdata <- function(object, newdata)
{
    # process newdata for models built with linmod.default

    x <- if(is.vector(newdata))     # allow newdata to be a vector
            matrix(newdata, ncol = length(object$coefficients) - 1)
         else
            as.matrix(newdata)      # allow newdata to be a data.frame

    cbind(1, x)                     # return data with an intercept column
}
process.newdata.formula <- function(object, newdata)
{
    # process newdata for models built with linmod.formula

    newdata <- as.data.frame(newdata) # allows newdata to be a matrix
    terms <- object$terms
    dataClasses <- attr(terms, "dataClasses")
    iresp <- attr(terms, "response")
    terms <- delete.response(terms)
    # na.action=na.pass because we will catch NAs after (for clearer error msg)
    # xlevels is needed to convert strings to factor levels, for example:
    #    a <- linmod(Sepal.Length~Species,data=iris)
    #    predict(a,newdata=data.frame(Species="setosa"))
    mf <- model.frame(terms, newdata, na.action = na.pass, xlev = object$xlevels)
    if(anyNA(mf))
        stop("NA in 'newdata'")
    if(NROW(mf) != NROW(newdata)) {
        # Get here when model.frame() issues
        #   Warning: 'newdata' had M rows but variables found have N rows
        # Must stop, else the call to model.matrix() below would silently return bad data.

        # If a variable is missing, print its name to help the user.
        # TODO This will erroneously identify "sqrt(x)" as a missing var in the
        #      formula "y ~ sqrt(x)" (because the var is wrapped in a func call).
        varnames <- names(dataClasses)
        varnames <- varnames[-iresp]
        missing <- which(!(varnames %in% colnames(newdata)))
        missing.msg <- ""
        if(length(missing))
            missing.msg <- paste0(" (variable '", varnames[missing[1]],
                                  "' may be missing from newdata)")
        stop("newdata has ", NROW(newdata),
             " rows but model.frame returned ", NROW(mf), " rows", missing.msg)
    }
    .checkMFClasses(dataClasses, mf) # check types in newdata match original data
    model.matrix(terms, mf)
}
do.predict.linmod <- function(object, x)
{
    # workhorse function for linear model predictions
    # processing by model.matrix etc. and sanity checks on x already completed
    # x is a numeric matrix (if model has intercept, first col of x is all 1s)

    x %*% coef(object)
}
summary.linmod <- function(object = stop("no 'object' argument"), ...)
{
    stop.if.dot.arg.used(...)
    se <- sqrt(diag(object$vcov))
    t.value <- coef(object) / se
    p.value <-
        if(object$df.residual == 0) # avoid warning from pt()
            rep_len(0, length.out=length(t.value))
        else
            2 * pt(-abs(t.value), df = object$df.residual)

    coefficients <- cbind(Estimate = coef(object),
                          StdErr   = se,
                          t.value  = t.value,
                          p.value  = p.value)

    retval <- list(call         = object$call,
                   coefficients = coefficients)

    class(retval) <- "summary.linmod"

    retval
}
print.linmod <- function(x = stop("no 'x' argument"), ...)
{
    stop.if.dot.arg.used(...)
    print.model.call(x)
    print(x$coefficients)
    invisible(x)
}
print.summary.linmod <- function(x = stop("no 'x' argument"), ...)
{
    stop.if.dot.arg.used(...)
    print.model.call(x)
    print(x$coefficients)
    invisible(x)
}
print.model.call <- function(x)
{
    cat("Call: ") # print.lm has a newline here, but a space is more compact
    # use paste0 to convert vector of strings to single string if necessary
    cat(strwrap(paste0(deparse(x$call, control = NULL, nlines = 5),
                       sep = " ", collapse = " "), exdent = 6), sep = "\n")
    cat("\n")
}
# stop.if.dot.arg.used will cause an error message if any args are passed to it.
# We use it to test if any dots arg of the calling function was used, for
# functions that must have a dots arg (to match the generic method) but don't
# actually use the dots.  This helps the user catch mistyped or illegal args.
# R version 3.3-0 or higher has a function chkDots which could be used instead.

stop.if.dot.arg.used <- function()
{
    NULL
}
