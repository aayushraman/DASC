#' Update G in Semi-NMF
#'
#' @param X Data expression matrix need to be factorized
#' @param mf The basis matrix
#' @param mg The co-efficient matrix
#' @return \code{G} The basis matrix
#' @author Haidong Yi, Ayush T. Raman
#' @export
#' @examples
#' X <- matrix(1:12,nrow=4)
#' mf <- matrix(1:8,nrow=4)
#' mg <- matrix(1:6,ncol=2)
#' mg <- updateG(X,mf,mg)

updateG <- function(X, mf, mg) {
    p <- nrow(X)
    n <- ncol(X)
    k <- ncol(mf)
    G <- matrix(nrow = n, ncol = k)
    GFFP <- t(mf) %*% mf
    GFFN <- t(mf) %*% mf
    XF <- t(X) %*% mf
    for (i in 1:k) {
        for (j in 1:k) {
            GFFP[i, j] <- 0.5 * (abs(GFFP[i, j]) + GFFP[i, j])
            GFFN[i, j] <- 0.5 * (abs(GFFN[i, j]) - GFFN[i, j])
        }
    }
    GFFP <- mg %*% GFFP
    GFFN <- mg %*% GFFN
    for (i in 1:n) {
        for (j in 1:k) {
            num <- 0.5 * (abs(XF[i, j]) + XF[i, j]) + GFFN[i, j]
            den <- 0.5 * (abs(XF[i, j]) - XF[i, j]) + GFFP[i, j]
            G[i, j] <- mg[i, j] * sqrt(num / den)
        }
    }
    G
}
