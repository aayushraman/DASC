#' Get Laplace matrix from factor vector
#'
#' @param colDataBATCH A factor vector
#' @return \code{adjacency} The Laplace matrix of \code{colDataBATCH}
#' @export
#' @author Haidong Yi, Ayush T. Raman
#' @examples
#' batch.factor = c(rep('human',13),rep('mouse',13))
#' batch.factor = as.factor(batch.factor)
#' adj = trans_Laplace(batch.factor)
#'

trans_Laplace <- function(colDataBATCH) {
    adjacency <- matrix(0, nrow = length(colDataBATCH), ncol = length(colDataBATCH))
    iter <- length(levels(colDataBATCH))
    for (i in 1:iter) {
        FACTOR <- which(colDataBATCH == levels(colDataBATCH)[i])
        N <- length(FACTOR)
        for (j in 1:(N - 1)) {
            for (k in (j + 1):N) {
                adjacency[FACTOR[j], FACTOR[k]] <- 1
            }
        }
    }
    adjacency <- adjacency + t(adjacency)
    adjacency <- diag(rowSums(adjacency)) - adjacency
    return(adjacency)
}
