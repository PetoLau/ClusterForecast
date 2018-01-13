source("kmeanspp2.R")

clusterOptimKmeans <- function(matrixOOM, k_low, k_high, parallel = F, ncores) {
  
  mat <- data.matrix(matrixOOM)

  if (parallel == FALSE) {
    clusterings <- lapply(c(k_low:k_high), function(x) kmeanspp(mat, x, nstart = 10, iter.max = 20))
  } else {
  # cl <- makeCluster(ncores)
  cl <- makeForkCluster(ncores, outfile = "")
  registerDoParallel(cl)
  # clusterExport(cl, varlist = c("mat", "kmeanspp"), envir = environment())
  clusterings <- parLapply(cl, c(k_low:k_high), function(x) kmeanspp(mat, x, nstart = 10, iter.max = 20))
  if(!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- c()
  }
  }
  
  if (parallel == FALSE) {
    DB_values <- sapply(seq_along(clusterings), function(x) intCriteria(mat, as.integer(clusterings[[x]]$cluster), c("Davies_Bouldin")))
  } else {
    # cl <- makeCluster(ncores)
    # clusterExport(cl, varlist = c("clusterings", "mat", "intCriteria"), envir = environment())
    cl <- makeForkCluster(ncores, outfile = "")
    registerDoParallel(cl)
    DB_values <- parSapply(cl, seq_along(clusterings), function(x) intCriteria(mat, as.integer(clusterings[[x]]$cluster), c("Davies_Bouldin")))
    if(!is.null(cl)){
      parallel::stopCluster(cl)
      cl <- c()
    }
  }
  
  return(list(clustering = clusterings[[which.min(DB_values)]]$cluster, centroids = clusterings[[which.min(DB_values)]]$centers))
}

clusterOptimKmedoids <- function(matrixOOM, k_low, k_high, parallel = F, ncores) {
  
  mat <- data.matrix(matrixOOM)
  
  if (parallel == FALSE) {
    clusterings <- lapply(c(k_low:k_high), function(x) pam(mat, x, metric = "euclidean"))
  } else {
    # cl <- makeCluster(ncores)
    # clusterExport(cl, varlist = c("mat", "pam"), envir = environment())
    cl <- makeForkCluster(ncores, outfile = "")
    registerDoParallel(cl)
    clusterings <- parLapply(cl, c(k_low:k_high), function(x) pam(mat, x, metric = "euclidean"))
    if(!is.null(cl)) {
      parallel::stopCluster(cl)
      cl <- c()
    }
  }

  if (parallel == FALSE) {
    DB_values <- sapply(seq_along(clusterings), function(x) intCriteria(mat, as.integer(clusterings[[x]]$clustering), c("Davies_Bouldin")))
  } else{
    # cl <- makeCluster(ncores)
    # clusterExport(cl, varlist = c("clusterings", "mat", "intCriteria"), envir = environment())
    cl <- makeForkCluster(ncores, outfile = "")
    registerDoParallel(cl)
    DB_values <- parSapply(cl, seq_along(clusterings), function(x) intCriteria(mat, as.integer(clusterings[[x]]$clustering), c("Davies_Bouldin")))
    if(!is.null(cl)){
      parallel::stopCluster(cl)
      cl <- c()
    }
  }

  return(list(clustering = clusterings[[which.min(DB_values)]]$clustering, centroids = clusterings[[which.min(DB_values)]]$medoids))
}

check.matrix <- function(mat) {
  if (is.matrix(mat) == TRUE) {
    return(mat)
  } else
    return(t(as.matrix(mat)))
}
