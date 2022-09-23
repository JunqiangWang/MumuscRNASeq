#' Nice Exponential representations of scientific notation
#' 
#' This function generates a plotmath or latex representation of scientific notation
#' 
#' @param x Numeric vector
#' @param drop.1 Logical, whether 1 in 1 x type of representatons should be dropped
#' @param sub10 Either logical, "10", a non-negative integer or a length 2 integer vector, indicating if some expression should be formatted traditionally, when integer, all expression before the integer are simplified. when a 2 elements vector, all between the indicated range are simplified
#' @param digits Number of significant digits
#' @param lab.type Character string indicating how the result should look like, either plotmath or latex
#' @param lab.sep Character separator between mantissa and exponent
#' @return Vector of formated numbers
#' @export
niceExponent <- function(x, drop.1 = TRUE, sub10 = "10", digits = 2, digits.fuzz, lab.type = c("plotmath", "latex"), lab.sep = c("cdot", "times"))
{
    lab.type <- match.arg(lab.type)
    lab.sep <- match.arg(lab.sep)    
    eT <- floor(log10(abs(x)) + 10^-digits)
    mT <- signif(x / 10^eT, digits)
    ss <- vector("list", length(x))
    if(sub.10 <- !identical(sub10, FALSE)) {
        if(identical(sub10, TRUE))
            sub10 <- c(0,0)
        else if(identical(sub10, "10"))
            sub10 <- 0:1
        sub10 <- as.integer(sub10)
        noE <-
            if(length(sub10) == 1) {
                if(sub10 < 0)
                    stop("'sub10' must not be negative if a single number")
                eT <= sub10
            } else if(length(sub10) == 2) {
                stopifnot(sub10[1] <= sub10[2])
                sub10[1] <= eT & eT <= sub10[2]
            } else stop("invalid 'sub10'")
        mT[noE] <- mT[noE] * 10^eT[noE]
    }
    if (lab.type == "plotmath") {
        for(i in seq(along = x))
            ss[[i]] <-
            if(x[i] == 0) quote(0)
        else if(sub.10 &&  noE[i]    ) substitute( A, list(A = mT[i]))
        else if(drop.1 && mT[i] ==  1) substitute( 10^E, list(E = eT[i]))
        else if(drop.1 && mT[i] == -1) substitute(-10^E, list(E = eT[i]))
        else substitute(A %*% 10^E, list(A = mT[i], E = eT[i]))
        do.call("expression", ss)
    }
    else { 
        mTf <- format(mT)
        eTf <- format(eT)
        for(i in seq(along = x))
            ss[[i]] <-
            if(x[i] == 0) ""
        else if(sub.10 &&  noE[i]    ) mTf[i]
        else if(drop.1 && mT[i] ==  1) sprintf("$10^{%s}$", eTf[i])
        else if(drop.1 && mT[i] == -1) sprintf("$-10^{%s}$",eTf[i])
        else sprintf("$%s \\%s 10^{%s}$", mTf[i], lab.sep,  eTf[i])
        ss  ## perhaps unlist(ss) ?
    }
}

