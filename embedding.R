reduceDimension <- function (cds, max_components = 2, 
                             reduction_method = c("DDRTree", "ICA", "tSNE", "SimplePPT", 
                                                  "L1-graph", "SGL-tree"), 
                             norm_method = c("log", "vstExprs", "none"), 
                             residualModelFormulaStr = NULL, pseudo_expr = 1, 
                             relative_expr = TRUE, auto_param_selection = TRUE, verbose = FALSE, 
                             scaling = TRUE, ...) {
  extra_arguments <- list(...)
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0, ]
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose) 
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), 
                                       data = pData(cds), drop.unused.levels = TRUE)
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  } else {
    X.model_mat <- NULL
  }
  if (scaling) {
    FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
    FM <- FM[!is.na(row.names(FM)), ]
  } else FM <- as.matrix(FM)
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ]
  if (is.function(reduction_method)) {
    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    reducedDimW(cds) <- as.matrix(reducedDim)
    reducedDimA(cds) <- as.matrix(reducedDim)
    reducedDimS(cds) <- as.matrix(reducedDim)
    reducedDimK(cds) <- as.matrix(reducedDim)
    dp <- as.matrix(dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  } else {
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "tSNE") {
      if (verbose) 
        message("Remove noise by PCA ...")
      if ("num_dim" %in% names(extra_arguments)) {
        num_dim <- extra_arguments$num_dim
      } else {
        num_dim <- 50
      }
      FM <- (FM)
      irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, 
                                               min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
      irlba_pca_res <- irlba_res$x
      topDim_pca <- irlba_pca_res
      if (verbose) 
        message("Reduce dimension by tSNE ...")
      tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, 
                               pca = F, ...)
      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)
      reducedDimA(cds) <- t(tsne_data)
      cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
      cds@auxClusteringData[["tSNE"]]$reduced_dimension <- t(topDim_pca)
      cds@dim_reduce_type <- "tSNE"
    } else if (reduction_method == "ICA") {
      if (verbose) 
        message("Reducing to independent components")
      init_ICA <- ica_helper(Matrix::t(FM), max_components, 
                             use_irlba = TRUE, ...)
      x_pca <- Matrix::t(Matrix::t(FM) %*% init_ICA$K)
      W <- Matrix::t(init_ICA$W)
      weights <- W
      A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
      colnames(A) <- colnames(weights)
      rownames(A) <- rownames(FM)
      S <- weights %*% x_pca
      rownames(S) <- colnames(weights)
      colnames(S) <- colnames(FM)
      reducedDimW(cds) <- as.matrix(W)
      reducedDimA(cds) <- as.matrix(A)
      reducedDimS(cds) <- as.matrix(S)
      reducedDimK(cds) <- as.matrix(init_ICA$K)
      adjusted_S <- Matrix::t(reducedDimS(cds))
      dp <- as.matrix(dist(adjusted_S))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "ICA"
    } else if (reduction_method == "DDRTree") {
      if (verbose) 
        message("Learning principal graph with DDRTree")
      if (auto_param_selection & ncol(cds) >= 100) {
        if ("ncenter" %in% names(extra_arguments)) 
          ncenter <- extra_arguments$ncenter
        else 
          ncenter <- cal_ncenter(ncol(FM))
        ddr_args <- c(list(X = FM, dimensions = max_components, 
                           ncenter = ncenter, verbose = verbose), extra_arguments[names(extra_arguments) %in% 
                                                                                    c("initial_method", "maxIter", "sigma", "lambda", 
                                                                                      "param.gamma", "tol")])
        ddrtree_res <- do.call(DDRTree, ddr_args)
      } else {
        ddrtree_res <- DDRTree(FM, max_components, verbose = verbose, 
                               ...)
      }
      if (ncol(ddrtree_res$Y) == ncol(cds)) 
        colnames(ddrtree_res$Y) <- colnames(FM)
      else colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      colnames(ddrtree_res$Z) <- colnames(FM)
      reducedDimW(cds) <- ddrtree_res$W
      reducedDimS(cds) <- ddrtree_res$Z
      reducedDimK(cds) <- ddrtree_res$Y
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- findNearestPointOnMST(cds)
    } else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}
DDRTree <- function (X, dimensions = 2, initial_method = NULL, maxIter = 20, 
                     sigma = 0.001, lambda = NULL, ncenter = NULL, param.gamma = 10, 
                     tol = 0.001, verbose = F, ...) {
  D <- nrow(X)
  N <- ncol(X)
  W <- pca_projection_R(X %*% t(X), dimensions)
  if (is.null(initial_method)) {
    Z <- t(W) %*% X
  } else {
    tmp <- initial_method(X, ...)
    if (ncol(tmp) > D | nrow(tmp) > N) 
      stop("The dimension reduction method passed need to return correct dimensions")
    Z <- tmp[, 1:dimensions]
    Z <- t(Z)
  }
  if (is.null(ncenter)) {
    K <- N
    Y <- Z[, 1:K]
  } else {
    K <- ncenter
    if (K > ncol(Z)) 
      stop("Error: ncenters must be greater than or equal to ncol(X)")
    centers = t(Z)[seq(1, ncol(Z), length.out = K), ]
    kmean_res <- kmeans(t(Z), K, centers = centers)
    Y <- kmean_res$centers
    Y <- t(Y)
  }
  if (is.null(lambda)) {
    lambda = 5 * ncol(X)
  }
  ddrtree_res <- DDRTree:::DDRTree_reduce_dim(X, Z, Y, W, dimensions, 
                                              maxIter, K, sigma, lambda, param.gamma, tol, verbose)
  return(list(W = ddrtree_res$W, Z = ddrtree_res$Z, stree = ddrtree_res$stree, 
              Y = ddrtree_res$Y, X = ddrtree_res$X, objective_vals = ddrtree_res$objective_vals, 
              history = NULL))
}



cal_ncenter <- function (ncells, ncells_limit = 100) {
  round(2 * ncells_limit * log(ncells)/(log(ncells) + log(ncells_limit)))
}

