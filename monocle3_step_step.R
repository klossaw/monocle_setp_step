cds
use_partition = TRUE
close_loop = TRUE
learn_graph_control = NULL
verbose = FALSE
# ====================== learn_graph ======================
{
  reduction_method <- "UMAP"
  if (!is.null(learn_graph_control)) {
    assertthat::assert_that(methods::is(learn_graph_control, 
                                        "list"))
    assertthat::assert_that(all(names(learn_graph_control) %in% 
                                  c("euclidean_distance_ratio", "geodesic_distance_ratio", 
                                    "minimal_branch_len", "orthogonal_proj_tip", 
                                    "prune_graph", "scale", "ncenter", "nn.k", "rann.k", 
                                    "maxiter", "eps", "L1.gamma", "L1.sigma", "nn.method", 
                                    "nn.metric", "nn.n_trees", "nn.search_k", "nn.M", 
                                    "nn.ef_construction", "nn.ef", "nn.grain_size", 
                                    "nn.cores")), msg = "Unknown variable in learn_graph_control")
  }
  if (!is.null(learn_graph_control[["rann.k"]]) && !is.null(learn_graph_control[["nn.k"]])) {
    assertthat::assert_that(learn_graph_control[["rann.k"]] == 
                              learn_graph_control[["nn.k"]], msg = paste0("both learn_graph_control$nn.k and learn_graph_control$rann.k are", 
                                                                          " defined and are unequal. See help(learn_graph) for more", 
                                                                          " information."))
  }
  if (is.null(learn_graph_control[["nn.k"]]) && !is.null(learn_graph_control[["rann.k"]])) 
    learn_graph_control[["nn.k"]] <- learn_graph_control[["rann.k"]]
  
  euclidean_distance_ratio <- ifelse(is.null(learn_graph_control$euclidean_distance_ratio), 
                                     1, learn_graph_control$euclidean_distance_ratio)
  geodesic_distance_ratio <- ifelse(is.null(learn_graph_control$geodesic_distance_ratio), 
                                    1/3, learn_graph_control$geodesic_distance_ratio)
  minimal_branch_len <- ifelse(is.null(learn_graph_control$minimal_branch_len), 
                               10, learn_graph_control$minimal_branch_len)
  orthogonal_proj_tip <- ifelse(is.null(learn_graph_control$orthogonal_proj_tip), 
                                FALSE, learn_graph_control$orthogonal_proj_tip)
  prune_graph <- ifelse(is.null(learn_graph_control$prune_graph), 
                        TRUE, learn_graph_control$prune_graph)
  ncenter <- learn_graph_control$ncenter
  scale <- ifelse(is.null(learn_graph_control$scale), FALSE, 
                  learn_graph_control$scale)
  nn.k <- ifelse(is.null(learn_graph_control[["nn.k"]]), 25, 
                 learn_graph_control[["nn.k"]])
  maxiter <- ifelse(is.null(learn_graph_control$maxiter), 10, 
                    learn_graph_control$maxiter)
  eps <- ifelse(is.null(learn_graph_control$eps), 1e-05, learn_graph_control$eps)
  L1.gamma <- ifelse(is.null(learn_graph_control$L1.gamma), 
                     0.5, learn_graph_control$L1.gamma)
  L1.sigma <- ifelse(is.null(learn_graph_control$L1.sigma), 
                     0.01, learn_graph_control$L1.sigma)
  
  
  
  nn_control <- list()
  if (!is.null(learn_graph_control[["nn.method"]])) 
    nn_control[["method"]] <- learn_graph_control[["nn.method"]]
  if (!is.null(learn_graph_control[["nn.metric"]])) 
    nn_control[["metric"]] <- learn_graph_control[["nn.metric"]]
  if (!is.null(learn_graph_control[["nn.n_trees"]])) 
    nn_control[["n_trees"]] <- learn_graph_control[["nn.n_trees"]]
  if (!is.null(learn_graph_control[["nn.search_k"]])) 
    nn_control[["search_k"]] <- learn_graph_control[["nn.search_k"]]
  if (!is.null(learn_graph_control[["nn.M"]])) 
    nn_control[["M"]] <- learn_graph_control[["nn.M"]]
  if (!is.null(learn_graph_control[["nn.ef_construction"]])) 
    nn_control[["ef_construction"]] <- learn_graph_control[["nn.ef_construction"]]
  if (!is.null(learn_graph_control[["nn.ef"]])) 
    nn_control[["ef"]] <- learn_graph_control[["nn.ef"]]
  if (!is.null(learn_graph_control[["nn.grain_size"]])) 
    nn_control[["grain_size"]] <- learn_graph_control[["nn.grain_size"]]
  if (!is.null(learn_graph_control[["nn.cores"]])) 
    nn_control[["cores"]] <- learn_graph_control[["nn.cores"]]
  if (verbose) 
    report_nn_control("nn_control: ", nn_control)
  
  
  nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
  nn_control <- set_nn_control(mode = 3, nn_control = nn_control, 
                               nn_control_default = nn_control_default, nn_index = NULL, 
                               k = nn.k, verbose = verbose)
  if (use_partition) {
    partition_list <- cds@clusters[[reduction_method]]$partitions
  }
  else {
    partition_list <- rep(1, nrow(colData(cds)))
  }
  multi_tree_DDRTree_res <- multi_component_RGE(cds, scale = scale, 
                                                reduction_method = reduction_method, partition_list = partition_list, 
                                                irlba_pca_res = SingleCellExperiment::reducedDims(cds)[[reduction_method]], 
                                                max_components = max_components, ncenter = ncenter, nn.k = nn.k, 
                                                nn_control = nn_control, maxiter = maxiter, eps = eps, 
                                                L1.gamma = L1.gamma, L1.sigma = L1.sigma, close_loop = close_loop, 
                                                euclidean_distance_ratio = euclidean_distance_ratio, 
                                                geodesic_distance_ratio = geodesic_distance_ratio, prune_graph = prune_graph, 
                                                minimal_branch_len = minimal_branch_len, verbose = verbose)
  rge_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
  rge_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
  rge_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
  cds <- multi_tree_DDRTree_res$cds
  dp_mst <- multi_tree_DDRTree_res$dp_mst
  principal_graph(cds)[[reduction_method]] <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- rge_res_Y
  cds <- project2MST(cds, project_point_to_line_segment, orthogonal_proj_tip, 
                     verbose, reduction_method, rge_res_Y)
  cds
}


# ====================== multi_component_RGE ======================

function (cds, scale = FALSE, reduction_method, partition_list, 
          max_components, ncenter, irlba_pca_res, nn.k = 25, nn_control = list(), 
          maxiter, eps, L1.gamma, L1.sigma, close_loop = FALSE, euclidean_distance_ratio = 1, 
          geodesic_distance_ratio = 1/3, prune_graph = TRUE, minimal_branch_len = minimal_branch_len, 
          verbose = FALSE) 
{
  cluster <- NULL
  X <- t(irlba_pca_res)
  dp_mst <- NULL
  pr_graph_cell_proj_closest_vertex <- NULL
  cell_name_vec <- NULL
  reducedDimK_coord <- NULL
  merge_rge_res <- NULL
  max_ncenter <- 0
  for (cur_comp in sort(unique(partition_list))) {
    if (verbose) {
      message("Processing partition component ", cur_comp)
    }
    X_subset <- X[, partition_list == cur_comp]
    if (verbose) 
      message("Current partition is ", cur_comp)
    if (scale) {
      X_subset <- t(as.matrix(scale(t(X_subset))))
    }
    if (is.null(ncenter)) {
      num_clusters_in_partition <- length(unique(clusters(cds, 
                                                          reduction_method)[colnames(X_subset)]))
      num_cells_in_partition = ncol(X_subset)
      curr_ncenter <- cal_ncenter(num_clusters_in_partition, 
                                  num_cells_in_partition)
      if (is.null(curr_ncenter) || curr_ncenter >= ncol(X_subset)) {
        curr_ncenter <- ncol(X_subset) - 1
      }
    }
    else {
      curr_ncenter <- min(ncol(X_subset) - 1, ncenter)
    }
    if (verbose) 
      message("Using ", curr_ncenter, " nodes for principal graph")
    kmean_res <- NULL
    centers <- t(X_subset)[seq(1, ncol(X_subset), length.out = curr_ncenter), 
                           , drop = FALSE]
    centers <- centers + matrix(stats::rnorm(length(centers), 
                                             sd = 1e-10), nrow = nrow(centers))
    kmean_res <- tryCatch({
      stats::kmeans(t(X_subset), centers = centers, iter.max = 100)
    }, error = function(err) {
      stats::kmeans(t(X_subset), centers = curr_ncenter, 
                    iter.max = 100)
    })
    if (kmean_res$ifault != 0) {
      message("kmeans returned ifault = ", kmean_res$ifault)
    }
    nearest_center <- find_nearest_vertex(t(kmean_res$centers), 
                                          X_subset, process_targets_in_blocks = TRUE)
    medioids <- X_subset[, unique(nearest_center)]
    reduced_dim_res <- t(medioids)
    mat <- t(X_subset)
    if (is.null(nn.k)) {
      k <- round(sqrt(nrow(mat))/2)
      k <- max(10, k)
    }
    else {
      k <- nn.k
    }
    if (verbose) 
      message("Finding kNN with ", k, " neighbors")
    nn_method <- nn_control[["method"]]
    dx <- search_nn_matrix(subject_matrix = mat, query_matrix = mat, 
                           k = min(k, nrow(mat) - 1), nn_control = nn_control, 
                           verbose = verbose)
    if (nn_method == "annoy" || nn_method == "hnsw") 
      dx <- swap_nn_row_index_point(nn_res = dx, verbose = verbose)
    nn.index <- dx$nn.idx[, -1]
    nn.dist <- dx$nn.dists[, -1]
    if (verbose) 
      message("Calculating the local density for each sample based on kNNs ...")
    rho <- exp(-rowMeans(nn.dist))
    mat_df <- as.data.frame(mat)
    tmp <- mat_df %>% tibble::rownames_to_column() %>% 
      dplyr::mutate(cluster = kmean_res$cluster, 
                    density = rho) %>% dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 1, wt = density) %>% dplyr::arrange(-dplyr::desc(cluster))
    medioids <- X_subset[, tmp$rowname]
    reduced_dim_res <- t(medioids)
    graph_args <- list(X = X_subset, C0 = medioids, maxiter = maxiter, 
                       eps = eps, L1.gamma = L1.gamma, L1.sigma = L1.sigma, 
                       verbose = verbose)
    rge_res <- do.call(calc_principal_graph, graph_args)
    names(rge_res)[c(2, 4, 5)] <- c("Y", "R", "objective_vals")
    stree <- rge_res$W
    stree_ori <- stree
    if (close_loop) {
      reduce_dims_old <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])[, 
                                                                                       partition_list == cur_comp]
      connect_tips_res <- connect_tips(cds, pd = colData(cds)[partition_list == 
                                                                cur_comp, ], R = rge_res$R, stree = stree, reducedDimK_old = rge_res$Y, 
                                       reducedDimS_old = reduce_dims_old, k = 25, nn_control = nn_control, 
                                       kmean_res = kmean_res, euclidean_distance_ratio = euclidean_distance_ratio, 
                                       geodesic_distance_ratio = geodesic_distance_ratio, 
                                       medioids = medioids, verbose = verbose)
      stree <- connect_tips_res$stree
    }
    if (prune_graph) {
      if (verbose) {
        message("Running graph pruning ...")
      }
      stree <- prune_tree(stree_ori, as.matrix(stree), 
                          minimal_branch_len = minimal_branch_len)
      rge_res$Y <- rge_res$Y[, match(row.names(stree), 
                                     row.names(stree_ori))]
      rge_res$R <- rge_res$R[, match(row.names(stree), 
                                     row.names(stree_ori))]
      medioids <- medioids[, row.names(stree)]
    }
    if (is.null(merge_rge_res)) {
      if (ncol(rge_res$Y) < 1) 
        warning("bad loop: ncol(rge_res$Y) < 1")
      colnames(rge_res$Y) <- paste0("Y_", 1:ncol(rge_res$Y))
      merge_rge_res <- rge_res
      colnames(merge_rge_res$X) <- colnames(X_subset)
      row.names(merge_rge_res$R) <- colnames(X_subset)
      if (ncol(merge_rge_res$Y) < 1) 
        warning("bad loop: ncol(merge_rge_res$Y) < 1")
      colnames(merge_rge_res$R) <- paste0("Y_", 1:ncol(merge_rge_res$Y))
      merge_rge_res$R <- list(merge_rge_res$R)
      merge_rge_res$stree <- list(stree)
      merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)
    }
    else {
      colnames(rge_res$X) <- colnames(X_subset)
      row.names(rge_res$R) <- colnames(X_subset)
      colnames(rge_res$R) <- paste0("Y_", (ncol(merge_rge_res$Y) + 
                                             1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), 
                                    sep = "")
      colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 
                                            1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), 
                                   sep = "")
      merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
      merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
      merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
      merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, 
                                        list(rge_res$objective_vals))
    }
    if (is.null(reducedDimK_coord)) {
      if (ncol(rge_res$Y) < 1) 
        warning("bad loop: ncol(rge_res$Y) < 1")
      curr_cell_names <- paste("Y_", 1:ncol(rge_res$Y), 
                               sep = "")
      pr_graph_cell_proj_closest_vertex <- matrix(apply(rge_res$R, 
                                                        1, which.max))
      cell_name_vec <- colnames(X_subset)
    }
    else {
      curr_cell_names <- paste("Y_", (ncol(reducedDimK_coord) + 
                                        1):(ncol(reducedDimK_coord) + ncol(rge_res$Y)), 
                               sep = "")
      pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, 
                                                 matrix(apply(rge_res$R, 1, which.max) + ncol(reducedDimK_coord)))
      cell_name_vec <- c(cell_name_vec, colnames(X_subset))
    }
    curr_reducedDimK_coord <- rge_res$Y
    dimnames(stree) <- list(curr_cell_names, curr_cell_names)
    cur_dp_mst <- igraph::graph.adjacency(stree, mode = "undirected", 
                                          weighted = TRUE)
    dp_mst <- igraph::graph.union(dp_mst, cur_dp_mst)
    reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
  }
  row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec
  ddrtree_res_W <- as.matrix(rge_res$W)
  ddrtree_res_Z <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  ddrtree_res_Y <- reducedDimK_coord
  R <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(cds), 
                                                          ncol(merge_rge_res$Y)))
  stree <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(merge_rge_res$Y), 
                                                              ncol(merge_rge_res$Y)))
  curr_row_id <- 1
  curr_col_id <- 1
  R_row_names <- NULL
  if (length(merge_rge_res$R) < 1) 
    warning("bad loop: length(merge_rge_res$R) < 1")
  for (i in 1:length(merge_rge_res$R)) {
    current_R <- merge_rge_res$R[[i]]
    stree[curr_col_id:(curr_col_id + ncol(current_R) - 1), 
          curr_col_id:(curr_col_id + ncol(current_R) - 1)] <- merge_rge_res$stree[[i]]
    curr_row_id <- curr_row_id + nrow(current_R)
    curr_col_id <- curr_col_id + ncol(current_R)
    R_row_names <- c(R_row_names, row.names(current_R))
  }
  row.names(R) <- R_row_names
  R <- R[colnames(cds), ]
  cds@principal_graph_aux[[reduction_method]] <- list(stree = stree, 
                                                      Q = merge_rge_res$Q, R = R, objective_vals = merge_rge_res$objective_vals, 
                                                      history = merge_rge_res$history)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- as.data.frame(pr_graph_cell_proj_closest_vertex)[colnames(cds), 
                                                                                                                                    , drop = FALSE]
  if (ncol(ddrtree_res_Y) < 1) 
    warning("bad loop: ncol(ddrtree_res_Y) < 1")
  colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), 
                                    sep = "")
  return(list(cds = cds, ddrtree_res_W = ddrtree_res_W, ddrtree_res_Z = ddrtree_res_Z, 
              ddrtree_res_Y = ddrtree_res_Y, dp_mst = dp_mst))
}
# ================= calc_principal_graph ==============
monocle3:::calc_principal_graph
function (X, C0, maxiter = 10, eps = 1e-05, L1.gamma = 0.5, L1.sigma = 0.01, 
          verbose = TRUE) 
{
  C <- C0
  K <- ncol(C)
  objs <- c()
  if (maxiter < 1) 
    warning("bad loop: maxiter < 1")
  for (iter in 1:maxiter) {
    norm_sq <- repmat(t(colSums(C^2)), K, 1)
    Phi <- norm_sq + t(norm_sq) - 2 * t(C) %*% C
    g <- igraph::graph.adjacency(Phi, mode = "lower", diag = TRUE, 
                                 weighted = TRUE)
    g_mst <- igraph::mst(g)
    stree <- igraph::get.adjacency(g_mst, attr = "weight", 
                                   type = "lower")
    stree_ori <- stree
    stree <- as.matrix(stree)
    stree <- stree + t(stree)
    W <- stree != 0
    obj_W <- sum(sum(stree))
    res = soft_assignment(X, C, L1.sigma)
    P <- res$P
    obj_P <- res$obj
    obj <- obj_W + L1.gamma * obj_P
    objs = c(objs, obj)
    if (verbose) 
      message("iter = ", iter, " obj = ", obj)
    if (iter > 1) {
      relative_diff = abs(objs[iter - 1] - obj)/abs(objs[iter - 
                                                           1])
      if (relative_diff < eps) {
        if (verbose) 
          message("eps = ", relative_diff, ", converge.")
        break
      }
      if (iter >= maxiter) {
        if (verbose) 
          message("eps = ", relative_diff, " reach maxiter.")
      }
    }
    C <- generate_centers(X, W, P, L1.gamma)
  }
  return(list(X = X, C = C, W = W, P = P, objs = objs))
}
# =============== connect_tips ================
function (cds, pd, R, stree, reducedDimK_old, reducedDimS_old, 
          k = 25, nn_control = nn_control, weight = FALSE, qval_thresh = 0.05, 
          kmean_res, euclidean_distance_ratio = 1, geodesic_distance_ratio = 1/3, 
          medioids, verbose = FALSE) 
{
  random_seed <- 0L
  reduction_method <- "UMAP"
  if (is.null(row.names(stree)) & is.null(row.names(stree))) {
    if (ncol(stree) < 1) 
      warning("bad loop: ncol(stree) < 1")
    dimnames(stree) <- list(paste0("Y_", 1:ncol(stree)), 
                            paste0("Y_", 1:ncol(stree)))
  }
  stree <- as.matrix(stree)
  stree[stree != 0] <- 1
  mst_g_old <- igraph::graph_from_adjacency_matrix(stree, mode = "undirected")
  if (is.null(kmean_res)) {
    tmp <- matrix(apply(R, 1, which.max))
    row.names(tmp) <- colnames(reducedDimS_old)
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
    data <- t(reducedDimS_old[, ])
    cluster_result <- louvain_clustering(data = data, pd = pd[, 
    ], weight = weight, nn_index = NULL, k = k, nn_control = nn_control, 
    louvain_iter = 1, random_seed = 0L, verbose = verbose)
    cluster_result$optim_res$membership <- tmp[, 1]
  }
  else {
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
    tip_pc_points_kmean_clusters <- sort(kmean_res$cluster[names(tip_pc_points)])
    data <- t(reducedDimS_old[, ])
    cluster_result <- louvain_clustering(data = data, pd = pd[row.names(data), 
    ], weight = weight, nn_index = NULL, k = k, nn_control = nn_control, 
    louvain_iter = 1, random_seed = random_seed, verbose = verbose)
    cluster_result$optim_res$membership <- kmean_res$cluster
  }
  cluster_graph_res <- compute_partitions(cluster_result$g, 
                                          cluster_result$optim_res, qval_thresh = qval_thresh, 
                                          verbose = verbose)
  dimnames(cluster_graph_res$cluster_mat) <- dimnames(cluster_graph_res$num_links)
  valid_connection <- which(cluster_graph_res$cluster_mat < 
                              qval_thresh, arr.ind = TRUE)
  valid_connection <- valid_connection[apply(valid_connection, 
                                             1, function(x) {
                                               all(x %in% tip_pc_points)
                                             }), ]
  G <- cluster_graph_res$cluster_mat
  G[cluster_graph_res$cluster_mat < qval_thresh] <- -1
  G[cluster_graph_res$cluster_mat > 0] <- 0
  G <- -G
  if (all(G == 0, na.rm = TRUE)) {
    return(list(stree = igraph::get.adjacency(mst_g_old), 
                Y = reducedDimK_old, G = G))
  }
  if (nrow(valid_connection) == 0) {
    return(list(stree = igraph::get.adjacency(mst_g_old), 
                Y = reducedDimK_old, G = G))
  }
  mst_g <- mst_g_old
  diameter_dis <- igraph::diameter(mst_g_old)
  reducedDimK_df <- reducedDimK_old
  pb4 <- utils::txtProgressBar(max = length(nrow(valid_connection)), 
                               file = "", style = 3, min = 0)
  res <- stats::dist(t(reducedDimK_old))
  g <- igraph::graph_from_adjacency_matrix(as.matrix(res), 
                                           weighted = TRUE, mode = "undirected")
  mst <- igraph::minimum.spanning.tree(g)
  max_node_dist <- max(igraph::E(mst)$weight)
  if (nrow(valid_connection) < 1) 
    warning("bad loop: nrow(valid_connection) < 1")
  for (i in 1:nrow(valid_connection)) {
    edge_vec <- sort(unique(cluster_result$optim_res$membership))[valid_connection[i, 
    ]]
    edge_vec_in_tip_pc_point <- igraph::V(mst_g_old)$name[edge_vec]
    if (length(edge_vec_in_tip_pc_point) == 1) 
      next
    if (all(edge_vec %in% tip_pc_points) & (igraph::distances(mst_g_old, 
                                                              edge_vec_in_tip_pc_point[1], edge_vec_in_tip_pc_point[2]) >= 
                                            geodesic_distance_ratio * diameter_dis) & (euclidean_distance_ratio * 
                                                                                       max_node_dist > stats::dist(t(reducedDimK_old[, edge_vec])))) {
      if (verbose) 
        message("edge_vec is ", edge_vec[1], "\t", edge_vec[2])
      if (verbose) 
        message("edge_vec_in_tip_pc_point is ", edge_vec_in_tip_pc_point[1], 
                "\t", edge_vec_in_tip_pc_point[2])
      mst_g <- igraph::add_edges(mst_g, edge_vec_in_tip_pc_point)
    }
    utils::setTxtProgressBar(pb = pb4, value = pb4$getVal() + 
                               1)
  }
  close(pb4)
  list(stree = igraph::get.adjacency(mst_g), Y = reducedDimK_df, 
       G = G)
}







# ====================== project2MST ======================
function (cds, Projection_Method, orthogonal_proj_tip = FALSE, 
          verbose, reduction_method, rge_res_Y) 
{
  target <- group <- distance_2_source <- rowname <- NULL
  dp_mst <- principal_graph(cds)[[reduction_method]]
  Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
  Y <- rge_res_Y
  cds <- findNearestPointOnMST(cds, reduction_method, rge_res_Y)
  closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex[, 1]]
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }
  else {
    if (length(Z) < 1) 
      warning("bad loop: length(Z) < 1")
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    if (length(Z[1:2, ]) < 1) 
      warning("bad loop: length(Z[1:2, ]) < 1")
    nearest_edges <- matrix(rep(0, length(Z[1:2, ])), ncol = 2)
    row.names(nearest_edges) <- colnames(cds)
    if (length(closest_vertex) < 1) 
      warning("bad loop: length(closest_vertex) < 1")
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(igraph::neighborhood(dp_mst, nodes = closest_vertex_names[i], 
                                              mode = "all")[[1]])[-1]
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          if (orthogonal_proj_tip) {
            tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                                              neighbor)])
          }
          else {
            tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                                                neighbor)])
          }
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                                              neighbor)])
        }
        if (any(is.na(tmp))) {
          tmp <- Y[, neighbor]
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, 
                                                  tmp)))
      }
      if (class(projection)[1] != "matrix") {
        projection <- as.matrix(projection)
      }
      which_min <- which.min(distance)
      P[, i] <- projection[which_min, ]
      nearest_edges[i, ] <- c(closest_vertex_names[i], 
                              neighbors[which_min])
    }
  }
  colnames(P) <- colnames(Z)
  dp_mst_list <- igraph::decompose.graph(dp_mst)
  dp_mst_df <- NULL
  partitions <- cds@clusters[[reduction_method]]$partitions
  assertthat::assert_that(!is.null(names(cds@clusters[[reduction_method]]$partitions)), 
                          msg = "names(cds@clusters[[reduction_method]]$partitions) == NULL")
  assertthat::assert_that(!(length(colnames(cds)) != length(names(cds@clusters[[reduction_method]]$partitions))), 
                          msg = "length( colnames(cds) ) != length(names(cds@clusters[[reduction_method]]$partitions))")
  assertthat::assert_that(!any(colnames(cds) != names(cds@clusters[[reduction_method]]$partitions)), 
                          msg = "colnames(cds)!=names(cds@clusters[[reduction_method]]$partitions)")
  if (length(dp_mst_list) == 1 & length(unique(partitions)) > 
      1) {
    partitions[partitions != "1"] <- "1"
  }
  if (!is.null(partitions)) {
    for (cur_partition in sort(unique(partitions))) {
      data_df <- NULL
      if (verbose) {
        message("\nProjecting cells to principal points for partition: ", 
                cur_partition)
      }
      subset_cds_col_names <- names(partitions[partitions == 
                                                 cur_partition])
      cur_z <- Z[, subset_cds_col_names]
      cur_p <- P[, subset_cds_col_names]
      if (ncol(cur_p) > 0 && nrow(cur_p) > 0) {
        cur_centroid_name <- igraph::V(dp_mst_list[[as.numeric(cur_partition)]])$name
        cur_nearest_edges <- nearest_edges[subset_cds_col_names, 
        ]
        data_df <- cbind(as.data.frame(t(cur_p)), apply(cur_nearest_edges, 
                                                        1, sort) %>% t())
        row.names(data_df) <- colnames(cur_p)
        if (nrow(cur_p) < 1) 
          warning("bad loop: nrow(cur_p) < 1")
        colnames(data_df) <- c(paste0("P_", 1:nrow(cur_p)), 
                               "source", "target")
        data_df$distance_2_source <- sqrt(colSums((cur_p - 
                                                     rge_res_Y[, as.character(data_df[, "source"])])^2))
        data_df <- data_df %>% tibble::rownames_to_column() %>% 
          dplyr::mutate(group = paste(source, target, 
                                      sep = "_")) %>% dplyr::arrange(group, dplyr::desc(-distance_2_source))
        data_df <- data_df %>% dplyr::group_by(group) %>% 
          dplyr::mutate(new_source = dplyr::lag(rowname), 
                        new_target = rowname)
        data_df[is.na(data_df$new_source), "new_source"] <- as.character(as.matrix(data_df[is.na(data_df$new_source), 
                                                                                           "source"]))
        added_rows <- which(is.na(data_df$new_source) & 
                              is.na(data_df$new_target))
        data_df <- as.data.frame(data_df, stringsAsFactors = FALSE)
        data_df <- as.data.frame(as.matrix(data_df), 
                                 stringsAsFactors = FALSE)
        data_df[added_rows, c("new_source", "new_target")] <- data_df[added_rows - 
                                                                        1, c("rowname", "target")]
        aug_P = cbind(cur_p, rge_res_Y)
        data_df$weight <- sqrt(colSums((aug_P[, data_df$new_source] - 
                                          aug_P[, data_df$new_target]))^2)
        data_df$weight <- data_df$weight + min(data_df$weight[data_df$weight > 
                                                                0])
        edge_list <- as.data.frame(igraph::get.edgelist(dp_mst_list[[as.numeric(cur_partition)]]), 
                                   stringsAsFactors = FALSE)
        dp <- as.matrix(stats::dist(t(rge_res_Y)[cur_centroid_name, 
        ]))
        edge_list$weight <- dp[cbind(edge_list[, 1], 
                                     edge_list[, 2])]
        colnames(edge_list) <- c("new_source", "new_target", 
                                 "weight")
        dp_mst_df <- Reduce(rbind, list(dp_mst_df, data_df[, 
                                                           c("new_source", "new_target", "weight")], edge_list))
      }
    }
  }
  dp_mst <- igraph::graph.data.frame(dp_mst_df, directed = FALSE)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_dist <- P
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}



