diff_test_helper <- function (x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, 
          relative_expr, weights, disp_func = NULL, verbose = FALSE) {
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, 
                                  sep = "")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, 
                               sep = "")
  
  
  x_orig <- x
  disp_guess <- 0
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (relative_expr == TRUE) {
      x <- x/Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE) {
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE) {
        if (expressionFamily@vfamily == "negbinomial") 
          expressionFamily <- negbinomial(isize = 1/disp_guess)
        else 
          expressionFamily <- negbinomial.size(size = 1/disp_guess)
      }
    }
  } else if (expressionFamily@vfamily %in% c("uninormal")) {
    f_expression <- x
  } else if (expressionFamily@vfamily %in% c("binomialff")) {
    f_expression <- x
  } else {
    f_expression <- log10(x)
  }
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")) {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon = 0.1, family = expressionFamily)
      } else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon = 0.1, family = expressionFamily))
      }
    } else {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon = 0.1, family = expressionFamily)
      } else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon = 0.1, family = expressionFamily))
      }
    }
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, error = function(e) {
    if (verbose) 
      print(e)
    data.frame(status = "FAIL", family = expressionFamily@vfamily, pval = 1, qval = 1)
  })
  test_res
}
smartEsApply <- function (X, MARGIN, FUN, convert_to_dense, ...) {
  parent <- environment(FUN)
  if (is.null(parent)) 
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  multiassign(names(pData(X)), pData(X), envir = e1)
  environment(FUN) <- e1
  if (monocle:::isSparseMatrix(exprs(X))) {
    res <- monocle:::sparseApply(exprs(X), MARGIN, FUN, convert_to_dense, 
                       ...)
  } else {
    res <- apply(exprs(X), MARGIN, FUN, ...)
  }
  if (MARGIN == 1) {
    names(res) <- row.names(X)
  } else {
    names(res) <- colnames(X)
  }
  res
}
sparseApply <- function (Sp_X, MARGIN, FUN, convert_to_dense, ...) {
  if (convert_to_dense) {
    if (MARGIN == 1) {
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[, i]), ...)
      }, FUN, ...)
    } else {
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[, i]), ...)
      }, FUN, ...)
    }
  } else {
    if (MARGIN == 1) {
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[, i], ...)
      }, FUN, ...)
    } else {
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[, i], ...)
      }, FUN, ...)
    }
  }
  return(res)
}
differentialGeneTest <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
          reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1, 
          verbose = FALSE) {
  status <- NA
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
  pd <- pData(cds)
  for (i in all_vars) {
    x <- pd[, i]
    if (any((c(Inf, NaN, NA) %in% x))) {
      stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
    }
  }
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  if (cores > 1) {
    diff_test_res <- mcesApply(cds, 1, diff_test_helper, 
                               c("BiocGenerics", "VGAM", "Matrix"), cores = cores, 
                               fullModelFormulaStr = fullModelFormulaStr, reducedModelFormulaStr = reducedModelFormulaStr, 
                               expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                               disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                               verbose = verbose)
    diff_test_res
  } else {
    diff_test_res <- smartEsApply(cds, 1, diff_test_helper, 
                                  convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr, 
                                  reducedModelFormulaStr = reducedModelFormulaStr, 
                                  expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                                  disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                                  verbose = verbose)
    diff_test_res
  }
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == "OK")] <- 
    p.adjust(subset(diff_test_res, status == "OK")[, "pval"], method = "BH")
  diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
  row.names(diff_test_res) <- diff_test_res[, 1]
  diff_test_res[, 1] <- NULL
  diff_test_res[row.names(cds), ]
}









BEAM <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
                  reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL, 
                  branch_point = 1, relative_expr = TRUE, branch_labels = NULL, 
                  verbose = FALSE, cores = 1, ...) {
  branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr, 
                               reducedModelFormulaStr = reducedModelFormulaStr, branch_states = branch_states, 
                               branch_point = branch_point, relative_expr = relative_expr, 
                               cores = cores, branch_labels = branch_labels, verbose = verbose, 
                               ...)
  cmbn_df <- branchTest_res[,1:4]
  if (verbose) 
    message("pass branchTest")
  fd <- fData(cds)[row.names(cmbn_df),]
  cmbn_df <- cbind(cmbn_df, fd)
  if (verbose) 
    message("return results")
  return(cmbn_df)
}
branchTest <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
                        reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL, 
                        branch_point = 1, relative_expr = TRUE, cores = 1, branch_labels = NULL, 
                        verbose = FALSE, ...) {
  if ("Branch" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
    cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = branch_states, 
                                         branch_point = branch_point, branch_labels = branch_labels, 
                                         ...)
  }
  else cds_subset <- cds
  branchTest_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, cores = cores, 
                                         relative_expr = relative_expr, verbose = verbose)
  return(branchTest_res)
}
buildBranchCellDataSet <- function (cds, progenitor_method = c("sequential_split", "duplicate"), 
                                    branch_states = NULL, branch_point = 1, branch_labels = NULL, 
                                    stretch = TRUE) {
  
  if (!is.null(branch_labels) & !is.null(branch_states)) {
    if (length(branch_labels) != length(branch_states)) 
      stop("length of branch_labels doesn't match with that of branch_states")
    branch_map <- setNames(branch_labels, as.character(branch_states))
  }
  if (cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  } else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }
  root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root_state <- pData(cds)[root_cell, ]$State
  pr_graph_root <- subset(pData(cds), State == root_state)
  if (cds@dim_reduce_type == "DDRTree") {
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
  } else {
    root_cell_point_in_Y <- row.names(pr_graph_root)
  }
  root_cell <- names(
    which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, mode = "all") == 1, useNames = T))[1]
  paths_to_root <- list()
  if (is.null(branch_states) == FALSE) {
    for (leaf_state in branch_states) {
      curr_cell <- subset(pData(cds), State == leaf_state)
      if (cds@dim_reduce_type == "DDRTree") {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell),]
      } else {
        curr_cell_point_in_Y <- row.names(curr_cell)
      }
      curr_cell <- names(which(degree(pr_graph_cell_proj_mst, 
                                      v = curr_cell_point_in_Y, 
                                      mode = "all") == 1, 
                               useNames = T))[1]
      path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, curr_cell, root_cell)
      path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
      if (cds@dim_reduce_type == "DDRTree") {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        ancestor_cells_for_branch <- row.names(closest_vertex)[
          which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_ancestor)]
      } else if (cds@dim_reduce_type == "ICA") {
        ancestor_cells_for_branch <- path_to_ancestor
      }
      ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, colnames(cds))
      paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
    }
  } else {
    if (cds@dim_reduce_type == "DDRTree") 
      pr_graph_cell_proj_mst <- minSpanningTree(cds)
    else 
      pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_cell <- mst_branch_nodes[branch_point]
    mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
    path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
    path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
    for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name) {
      descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], unreachable = FALSE)
      descendents <- descendents$order[!is.na(descendents$order)]
      descendents <- V(mst_no_branch_point)[descendents]$name
      if (root_cell %in% descendents == FALSE) {
        path_to_root <- unique(c(path_to_ancestor, branch_cell, descendents))
        if (cds@dim_reduce_type == "DDRTree") {
          closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
          path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
        } else {
          path_to_root <- path_to_root
        }
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        path_to_root <- intersect(path_to_root, colnames(cds))
        paths_to_root[[backbone_nei]] <- path_to_root
      }
    }
  }
  all_cells_in_subset <- c()
  if (is.null(branch_labels) == FALSE) {
    if (length(branch_labels) != 2) 
      stop("Error: branch_labels must have exactly two entries")
    names(paths_to_root) <- branch_labels
  }
  for (path_to_ancestor in paths_to_root) {
    if (length(path_to_ancestor) == 0) {
      stop("Error: common ancestors between selected State values on path to root State")
    }
    all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
  }
  all_cells_in_subset <- unique(all_cells_in_subset)
  common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
  cds <- cds[, row.names(pData(cds[, all_cells_in_subset]))]
  Pseudotime <- pData(cds)$Pseudotime
  pData <- pData(cds)
  if (stretch) {
    max_pseudotime <- -1
    for (path_to_ancestor in paths_to_root) {
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)
      if (max_pseudotime < max_pseudotime_on_path) {
        max_pseudotime <- max_pseudotime_on_path
      }
    }
    branch_pseudotime <- max(pData[common_ancestor_cells,]$Pseudotime)
    for (path_to_ancestor in paths_to_root) {
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)
      path_scaling_factor <- (max_pseudotime - branch_pseudotime)/(max_pseudotime_on_path - branch_pseudotime)
      if (is.finite(path_scaling_factor)) {
        branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
        pData[branch_cells,]$Pseudotime <- ((pData[branch_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
      }
    }
    pData$Pseudotime <- 100 * pData$Pseudotime/max_pseudotime
  }
  pData$original_cell_id <- row.names(pData)
  pData$original_cell_id <- row.names(pData)
  if (length(paths_to_root) != 2) 
    stop("more than 2 branch states are used!")
  pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1]
  progenitor_pseudotime_order <- order(pData[common_ancestor_cells, "Pseudotime"])
  if (progenitor_method == "duplicate") {
    ancestor_exprs <- exprs(cds)[, common_ancestor_cells]
    expr_blocks <- list()
    for (i in 1:length(paths_to_root)) {
      if (nrow(ancestor_exprs) == 1) 
        exprs_data <- t(as.matrix(ancestor_exprs))
      else exprs_data <- ancestor_exprs
      colnames(exprs_data) <- paste("duplicate", i, 1:length(common_ancestor_cells), sep = "_")
      expr_lineage_data <- exprs(cds)[, setdiff(paths_to_root[[i]], common_ancestor_cells)]
      exprs_data <- cbind(exprs_data, expr_lineage_data)
      expr_blocks[[i]] <- exprs_data
    }
    ancestor_pData_block <- pData[common_ancestor_cells,]
    pData_blocks <- list()
    weight_vec <- c()
    for (i in 1:length(paths_to_root)) {
      weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
      weight_vec_block <- rep(1, length(common_ancestor_cells))
      new_pData_block <- ancestor_pData_block
      row.names(new_pData_block) <- paste("duplicate", i, 1:length(common_ancestor_cells), sep = "_")
      pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], common_ancestor_cells),]
      weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
      weight_vec <- c(weight_vec, weight_vec_block)
      new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
      new_pData_block$Branch <- names(paths_to_root)[i]
      pData_blocks[[i]] <- new_pData_block
    }
    pData <- do.call(rbind, pData_blocks)
    exprs_data <- do.call(cbind, expr_blocks)
  } else if (progenitor_method == "sequential_split") {
    pData$Branch <- names(paths_to_root)[1]
    branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchA], "Branch"] <- names(paths_to_root)[1]
    branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchB], "Branch"] <- names(paths_to_root)[2]
    zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
    exprs_data <- cbind(exprs(cds), duplicate_root = exprs(cds)[, zero_pseudotime_root_cell])
    pData <- rbind(pData, pData[zero_pseudotime_root_cell,])
    row.names(pData)[nrow(pData)] <- "duplicate_root"
    pData[nrow(pData), "Branch"] <- names(paths_to_root)[2]
    weight_vec <- rep(1, nrow(pData))
    for (i in 1:length(paths_to_root)) {
      path_to_ancestor <- paths_to_root[[i]]
      branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
      pData[branch_cells, ]$Branch <- names(paths_to_root)[i]
    }
  }
  pData$Branch <- as.factor(pData$Branch)
  pData$State <- factor(pData$State)
  Size_Factor <- pData$Size_Factor
  fData <- fData(cds)
  colnames(exprs_data) <- row.names(pData)
  cds_subset <- newCellDataSet(as.matrix(exprs_data), 
                               phenoData = new("AnnotatedDataFrame", data = pData), 
                               featureData = new("AnnotatedDataFrame", data = fData), 
                               expressionFamily = cds@expressionFamily, 
                               lowerDetectionLimit = cds@lowerDetectionLimit)
  pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
  pData(cds_subset)$Size_Factor <- Size_Factor
  cds_subset@dispFitInfo <- cds@dispFitInfo
  return(cds_subset)
}





project2MST <- function (cds, Projection_Method) {
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(igraph::V(dp_mst)[suppressWarnings(nei(closest_vertex_names[i], mode = "all"))])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- monocle:::projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      if ("matrix" %in% class(projection)) 
        projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], 
      ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- as.matrix(dist(t(P)))
  min_dist = min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  cellPairwiseDistances(cds) <- dp
  gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- igraph::minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}
orderCells <- function (cds, root_state = NULL, num_paths = NULL, reverse = NULL) {
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  root_cell <- select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths = 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- as.matrix(dist(adjusted_S))
    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- extract_good_branched_ordering(res$subtree, 
                                                 res$root, 
                                                 cellPairwiseDistances(cds), num_paths, 
                                                 FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    minSpanningTree(cds) <- as.undirected(order_list$cell_ordering_tree)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)),]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
    minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- as.undirected(order_list$cell_ordering_tree)
  }
  else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- monocle:::extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),]$pseudo_time
    K_old <- reducedDimK(cds)
    old_dp <- cellPairwiseDistances(cds)
    old_mst <- minSpanningTree(cds)
    old_A <- reducedDimA(cds)
    old_W <- reducedDimW(cds)
    cds <- project2MST(cds, project_point_to_line_segment)
    minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- monocle:::extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)),]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(cds)$State <- cc_ordering[closest_vertex[,1], ]$cell_state
    }
    reducedDimK(cds) <- K_old
    cellPairwiseDistances(cds) <- old_dp
    minSpanningTree(cds) <- old_mst
    reducedDimA(cds) <- old_A
    reducedDimW(cds) <- old_W
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  }
  else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- monocle:::extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)),]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}


select_root_cell <- function (cds, root_state = NULL, reverse = FALSE) {
  if (is.null(root_state) == FALSE) {
    if (is.null(pData(cds)$State)) {
      stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
    }
    root_cell_candidates <- subset(pData(cds), State == root_state)
    if (nrow(root_cell_candidates) == 0) {
      stop(paste("Error: no cells for State =", root_state))
    }
    dp <- as.matrix(dist(t(reducedDimS(cds)[, row.names(root_cell_candidates)])))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 
                                1))
    diameter <- get.diameter(dp_mst)
    if (length(diameter) == 0) {
      stop(paste("Error: no valid root cells for State =", 
                 root_state))
    }
    root_cell_candidates <- root_cell_candidates[names(diameter), 
    ]
    if (is.null(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) == 
        FALSE && pData(cds)[cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell, 
        ]$State == root_state) {
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == 
                                                           min(root_cell_candidates$Pseudotime))]
    }
    else {
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == 
                                                           max(root_cell_candidates$Pseudotime))]
    }
    if (length(root_cell) > 1) 
      root_cell <- root_cell[1]
    if (cds@dim_reduce_type == "DDRTree") {
      graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell, 
      ]
      root_cell = V(minSpanningTree(cds))[graph_point_for_root_cell]$name
    }
  }
  else {
    if (is.null(minSpanningTree(cds))) {
      stop("Error: no spanning tree found for CellDataSet object. Please call reduceDimension before calling orderCells()")
    }
    diameter <- get.diameter(minSpanningTree(cds))
    if (is.null(reverse) == FALSE && reverse == TRUE) {
      root_cell = names(diameter[length(diameter)])
    }
    else {
      root_cell = names(diameter[1])
    }
  }
  return(root_cell)
}

extract_ddrtree_ordering <- function (cds, root_cell, verbose = T) {
  dp <- cellPairwiseDistances(cds)
  dp_mst <- minSpanningTree(cds)
  curr_state <- 1
  res <- list(subtree = dp_mst, root = root_cell)
  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name
  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name
  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name
  mst_traversal <- graph.dfs(dp_mst, root = root_cell, neimode = "all", 
                             unreachable = FALSE, father = TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  for (i in 1:length(mst_traversal$order)) {
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name
    if (is.na(mst_traversal$father[curr_node]) == FALSE) {
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, 
                                                         parent_node_name]
      if (degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
    }
    else {
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }
    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  ordering_df <- data.frame(sample_name = names(states), cell_state = factor(states), 
                            pseudo_time = as.vector(pseudotimes), parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

projPointOnLine <- function (point, line) {
  ap <- point - line[, 1]
  ab <- line[, 2] - line[, 1]
  res <- line[, 1] + c((ap %*% ab)/(ab %*% ab)) * ab
  return(res)
}
Projection_Method <- function (p, df) {
  A <- df[, 1]
  B <- df[, 2]
  AB <- (B - A)
  AB_squared = sum(AB^2)
  if (AB_squared == 0) {
    q <- A
  } else {
    Ap <- (p - A)
    t <- sum(Ap * AB)/AB_squared
    if (t < 0) {
      q <- A
    } else if (t > 1) {
      q <- B
    } else {
      q <- A + t * AB
    }
  }
  return(q)
}































