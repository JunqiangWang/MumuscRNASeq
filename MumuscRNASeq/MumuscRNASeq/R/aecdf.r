#' Approximate empirical commulative distribution function
#'
#' This function generates an empirical null model that computes a normalized statistics and p-value
#'
#' @param dnull Numerical vector representing the null model
#' @param symmetric Logical, whether the distribution should betreated as symmetric around zero and only one tail should be approximated
#' @param n Integer indicating the number of points to evaluate the empirical cummulative probability function
#' @return function with two parameters, \code{x} and \code{alternative}
#' @export
aecdf <- function(dnull, symmetric=FALSE, n=100) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        tmp <- sort(abs(dnull), decreasing=T)
        i <- 4
        n <- 4
        while(n<14) {
            i <- i+1
            n <- length(unique(tmp[1:i]))
            if (n==5) iq1 <- i
        }
        tl1 <- i
        iqr <- quantile(abs(dnull), c(.5, 1-iq1/length(dnull)))
        epd <- ecdf(abs(dnull))
        a <- list(x=knots(epd), y=epd(knots(epd)))
        fit <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl1:iq1)+1]-iqr[2], y=log(1-epd(iqr[2]))-log(1-a$y[length(a$x)-(tl1:iq1)+1])))
        val <- seq(0, iqr[2], length=n)
        pd <- approxfun(val, epd(val), method="linear", yleft=0, rule=2)
        dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
            alternative <- match.arg(alternative)
            x1 <- abs(x)
            p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
            p[!is.finite(p)] <- 1
            p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
            nes <- qnorm(p/2, lower.tail=F)*sign(x)
            switch(alternative,
                   two.sided={p <- p},
                   greater={p <- p/2; p[x<0] <- 1-p[x<0]},
                   less={p <- p/2; p[x>0] <- 1-p[x>0]}
            )
            names(nes) <- names(p) <- names(x)
            list(nes=nes, p.value=p)
        }
        return(dnull)
    }
    tmp <- sort(dnull, decreasing=FALSE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq1 <- i
    }
    tl1 <- i
    tmp <- sort(dnull, decreasing=TRUE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq2 <- i
    }
    tl2 <- i
    iqr <- quantile(dnull, c(iq1/length(dnull), .5, 1-iq2/length(dnull)))
    epd <- ecdf(dnull)
    a <- list(x=knots(epd), y=epd(knots(epd)))
    fit1 <- lm(y~0+x, data=list(x=a$x[iq1:tl1]-iqr[1], y=log(epd(iqr[1]))-log(a$y[iq1:tl1])))
    fit2 <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl2:iq2)+1]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(tl2:iq2)+1])))
    val <- seq(iqr[1], iqr[3], length=n)
    pd <- approxfun(val, epd(val), method="linear", rule=2)
    dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
        alternative <- match.arg(alternative)
        p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
        p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
        p1[!is.finite(p1)] <- 1
        p2[!is.finite(p2)] <- 1
        p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
        nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
        switch(alternative,
               two.sided={p <- p*2},
               greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
               less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
        )
        names(nes) <- names(p) <- names(x)
        list(nes=nes, p.value=p)
    }
    return(dnull)
}

#' Integration based on CDF (AOC)
#'
#' This function integrates a distribution of scores based on the area over the CDF curve
#'
#' @param x Numeric vector, matrix or list of vectors or matrixes
#' @param xlim Numeric vector of 2 elements indicating the range where to perform the integration
#' @details This function computes the area over the curve for the vector o columns of the matrix provided as input
#' @export

cdfInteg <- function(x, xlim=NULL) {
    if (is.null(xlim)) xlim <- range(unlist(x, use.names=FALSE), na.rm=TRUE)
    if (is.list(x)) return(sapply(x, cdfInteg, xlim=xlim))
    if (is.matrix(x)) return(apply(x, 2, cdfInteg, xlim=xlim))
    1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim)
}
