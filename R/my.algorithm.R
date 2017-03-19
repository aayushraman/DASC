#' Main function of semi-NMF
#'
#' @import NMF
#' @param target gene expression matrix
#' @param model Object of class: NMFfit
#' @param iternum Number of iterations
#' @return \code{model} the result from the semi-NMF algorithm
#' @author Haidong Yi, Ayush T. Raman

my.algorithm <- function(target, model, iternum = 100) {
    n <- 0
    while (TRUE) {
        ans <- t(updateG(target, basis(model), t(coef(model))))
        if (length(which(ans < 1e-20) > 0))
            break
        coef(model) <- ans
        basis(model) <- target %*% t(coef(model)) %*% 
                            solve(coef(model) %*% t(coef(model)))
        if (n > iternum || my_objective_function(model, target) < .5)
            break
        n <- n + 1
    }
    model
}
