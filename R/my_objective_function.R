#' Get the error of Semi-NMF
#'
#' @details
#' This is a customerized function defined in terms of \code{\link[NMF]{nmf}}. 
#' For more information, please go through the NMF vignette 
#' \url{https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf}
#'
#' @import NMF
#' @param model Object of class: NMFfit
#' @param target gene expression matrix
#' @return The result of semi-NMF for the current iteration
#'
#' @author Haidong Yi, Ayush T. Raman

my_objective_function <- function(model, target) {
    H <- NMF::basis(model)
    G <- t(coef(model))
    norm((target - H %*% t(G)), type = "f")
}
