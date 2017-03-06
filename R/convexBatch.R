#' Transform the adjacency matrix to a vector
#'
#' @param Adjacency the adjacency matrix of factor
#' @param n number of samples
#' @return \code{w} the vector of the adjacency matrix
#' @author Haidong Yi, Ayush T. Raman
#' @export
#' @examples
#' W = matrix(c(0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0),nrow=4)
#' w = W_toVector(W,4)

W_toVector <- function(Adjacency, n) {
    w <- double(n * (n - 1)/2)
    iter <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            w[iter] <- Adjacency[i, j]
            iter <- iter + 1
        }
    }
    return(w)
}

#' Representing node in this subtype
#'
#' @param v the index of the node
#' @param X the saved vector with the information of the parent of every node
#' @return \code{r} the parent index of the node
#' @author Haidong Yi, Ayush T. Raman

getFather <- function(v, X) {
    r <- v
    while (X[r] != r) {
        r <- X[r]
    }
    i <- v
    j <- numeric()
    while (i != r) {
        j <- X[i]
        X[i] <- r
        i <- j
    }
    return(r)
}

#' Merge two nodes
#'
#' @param x the index of the node
#' @param y the index of the node
#' @param X the saved vector with the information of the parent of every node
#' @return \code{X} A updated X vector with updates on father of every node
#' @author Haidong Yi, Ayush T. Raman

merge <- function(x, y, X) {
    fx <- getFather(x, X)
    fy <- getFather(y, X)
    if (fx < fy) {
        X[fx] <- fy
    } else {
        X[fy] <- fx
    }
    return(X)
}


#' Spanning tree from adjacency matrix
#'
#' @param ADJ the adjacency matrix of the factor
#' @return \code{ADJ} the spaning tree of the adjacency matrix
#' @author Haidong Yi, Ayush T. Raman
#' @export
#' @examples
#' W = matrix(c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0),nrow=4)
#' getclass(W)

getclass <- function(ADJ) {
    Rownum <- nrow(ADJ)
    Colnum <- ncol(ADJ)
    
    father <- numeric()
    for (i in 1:Rownum) {
        father[i] <- i
    }
    
    for (i in 2:Rownum) {
        for (j in 1:(i - 1)) {
            if (ADJ[i, j] > 0) {
                father <- merge(i, j, father)
            }
        }
    }
    sum <- numeric()
    pre <- numeric()
    j <- 1
    for (i in 1:Rownum) {
        pre[i] <- getFather(i, father)
        if (pre[i] == i) {
            sum[j] <- i
            j <- j + 1
        }
    }
    ADJ <- matrix(0, Rownum, Colnum)
    for (i in 1:Colnum) {
        if (i != pre[i]) {
            ADJ[i, pre[i]] <- 1
            ADJ[pre[i], i] <- 1
        }
    }
    ADJ
}

#' Batch detection via Semi-NMF
#'
#' @param edata the normalized target matrix, a data.frame The row is gene, the column is sample
#' @param pdata Phenotypic data summarizes information about the samples
#' @param factor A factor vector which controls the convex clustering
#' @param method Algorithm to use: 'admm' or 'ama'
#' @param type An integer indicating the norm used: 1 = 1-norm 2 = 2-norm 3 = 2-norm^2
#' @param lambda A double number The regularization parameter in the convex optimization
#' @param rank A integer sequence
#' @param nrun the iteration numbers of Semi-NMF
#' @param spanning parameter is assigned as false
#' @param annotation An annotation of the dataset
#' @return outputs the result of semi-NMF. It classifies each sample to its batch factor.
#' @import Biobase
#' @import cvxclustr
#' @import NMF
#' @import DESeq2
#' @import ggplot2
#' @import pcaExplorer
#' @export
#' @author Haidong Yi, Ayush T. Raman
#' @seealso \code{\link[cvxclustr]{cvxclust_path_ama}} and \code{\link[cvxclustr]{cvxclust_path_admm}} for
#' the detailed algorithm
#'

convexBatch <- function(edata, pdata, factor, method = "ama", type = 3, lambda, rank, nrun, spanning = FALSE, 
    annotation) {
    if (!is.null(type) && !(type %in% c(1, 2, 3))) 
        stop("type must be '1', '2', '3', or NULL")
    if (!is.null(method) && !(method %in% c("ama", "admm"))) 
        stop("method must be 'ama', 'admm', or NULL")
    edata = as.matrix(edata)
    
    Zero <- apply(edata, 1, sd)  #remove the gene without variance
    Zero.num <- which(Zero < 0.001)
    if (length(Zero.num) > 0) {
        names(Zero.num) <- NULL
        edata <- edata[-Zero.num, ]
    }
    edata = log(1 + edata)
    
    if (type == 3) {
        Laplace <- trans_Laplace(as.factor(factor))
        Udata <- edata %*% solve(diag(ncol(edata)) + lambda * Laplace)
        Udata = as.matrix(Udata)
        Bdata = edata - Udata
    } else {
        ADJ = trans_ADJ(as.factor(factor))
        if (spanning) {
            ADJ = getclass(ADJ)
        }
        w = W_toVector(ADJ, nrow(ADJ))
        sol = cvxclust(edata, w, lambda, method = method, type = type)
        Bdata = edata - sol$U[[1]]
    }
    
    Zero <- apply(Bdata, 1, sd)  #remove the gene without variance
    Zero.num <- which(Zero < 0.001)
    if (length(Zero.num) > 0) {
        names(Zero.num) <- NULL
        Bdata <- Bdata[-Zero.num, ]
    }
    
    metadata <- data.frame(labelDescription = names(pdata), row.names = names(pdata))
    pdata <- new("AnnotatedDataFrame", data = pdata, varMetadata = metadata)
    
    # get the ExpressionSet Class data for NMF calculate
    DataSet <- ExpressionSet(assayData = edata, phenoData = pdata, annotation = annotation)
    
    # Remark: In windows, the nrun > 1 will be crashed, So I recommend to use mac os x or linux to run the
    # code in R
    data.nmf <- nmf(DataSet, rank, my.algorithm, nrun = nrun, .opt = "v", objective = my_objective_function, 
        seed = my.seeding.method, mixed = TRUE)
    return(data.nmf)
}
