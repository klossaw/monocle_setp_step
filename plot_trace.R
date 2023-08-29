plot_cell_trajectory <- function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
                                  show_backbone = TRUE, backbone_color = "black", markers = NULL, 
                                  use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
                                  show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
                                  cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
                                  theta = 0, ...) {
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(cds)
  }
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
    select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>% 
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
    select_(source = "from", target = "to") %>% 
    left_join(
      ica_space_df %>% 
        select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
    left_join(
      ica_space_df %>% 
        select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                target_prin_graph_dim_2 = "prin_graph_dim_2"), 
      by = "target")
  data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
    select_(data_dim_1 = x, data_dim_2 = y) %>% rownames_to_column("sample_name") %>% 
    mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
                                         rownames_to_column("sample_name"), by = "sample_name")
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", by.y = "cell_id")
    if (use_color_gradient) {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
          geom_point(aes(color = value), size = I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~ feature_label)
      } else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
          geom_point(aes(color = log10(value + 0.1)), size = I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label)
      }
    } else {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, size = (value * 0.1))) + 
          facet_wrap(~feature_label)
      } else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, size = log10(value + 0.1))) + 
          facet_wrap(~feature_label)
      }
    }
  } else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    if (use_color_gradient) {
    } else {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    }
  } else {
    if (use_color_gradient) {
    } else {
      g <- g + geom_point(aes_string(color = color_by), size = I(cell_size), na.rm = TRUE)
    }
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>% 
      slice(match(mst_branch_nodes, sample_name)) %>% 
      mutate(branch_point_idx = seq_len(n()))
    g <- g + geom_point(
      aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), 
      size = 5, na.rm = TRUE, branch_point_df) + 
      geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", 
                           label = "branch_point_idx"), size = 4, 
                color = "white", na.rm = TRUE, branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_state_number) {
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  g <- g + monocle_theme_opts() + xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) + 
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) + 
    theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  g
}
plot_complex_cell_trajectory <- function (cds, x = 1, y = 2, root_states = NULL, color_by = "State", 
          show_tree = TRUE, show_backbone = TRUE, backbone_color = "black", 
          markers = NULL, show_cell_names = FALSE, cell_size = 1.5, 
          cell_link_size = 0.75, cell_name_size = 2, show_branch_points = TRUE, 
          ...) {
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", "SGL-tree")) {
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  if (is.null(reduced_dim_coords)) {
    stop("You must first call reduceDimension() before using this function")
  }
  dp_mst <- minSpanningTree(cds)
  if (is.null(root_states)) {
    if (is.null(lib_info_with_pseudo$Pseudotime)) {
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 1][1]
    } else root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))
    if (cds@dim_reduce_type != "ICA") 
      root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell,]]
  } else {
    candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
    if (cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, candidate_root_cells) == 1)]
    } else {
      Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells,]]
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, Y_candidate_root_cells) == 1)]
    }
  }
  tree_coords <- layout_as_tree(dp_mst, root = root_cell)
  ica_space_df <- data.frame(tree_coords)
  row.names(ica_space_df) <- colnames(reduced_dim_coords)
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
  ica_space_df$sample_name <- row.names(ica_space_df)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", by.y = "source", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", prin_graph_dim_2 = "source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, 
                   ica_space_df[, c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], 
                   by.x = "target", by.y = "sample_name", all = TRUE)
  edge_df <- plyr::rename(edge_df, 
                          c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                            prin_graph_dim_2 = "target_prin_graph_dim_2"))
  if (cds@dim_reduce_type == "ICA") {
    S_matrix <- tree_coords[, ]
  } else if (cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", "SGL-tree")) {
    S_matrix <- tree_coords[closest_vertex, ]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "row.names")
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", by.y = "cell_id")
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, I(cell_size))) + facet_wrap(~feature_label)
  } else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    if (class(data_df[, color_by]) == "numeric") {
      g <- g + geom_jitter(aes_string(
        color = paste0("log10(", color_by, " + 0.1)")), size = I(cell_size), na.rm = TRUE, 
        height = 5) + scale_color_viridis(name = paste0("log10(", color_by, ")"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), size = I(cell_size), na.rm = TRUE, height = 5)
    }
  } else {
    if (class(data_df[, color_by]) == "numeric") {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), 
                           size = I(cell_size), na.rm = TRUE, 
                           height = 5) + 
        scale_color_viridis(name = paste0("log10(", color_by, " + 0.1)"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), 
                           size = I(cell_size), na.rm = TRUE, height = 5)
    }
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx),]
    g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2"), 
                        size = 2 * cell_size, 
                        na.rm = TRUE, data = branch_point_df) + 
      geom_text(aes_string(x = "source_prin_graph_dim_1", 
                           y = "source_prin_graph_dim_2", 
                           label = "branch_point_idx"), 
                size = 1.5 * cell_size, color = "white", na.rm = TRUE, 
                data = branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + theme(strip.background = element_rect(colour = "white",fill = "white")) + 
    theme(panel.border = element_blank()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(legend.key = element_blank()) + xlab("") + ylab("") + 
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) + 
    theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(line = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  g
}
plot_multiple_branches_pseudotime <- function (cds, branches, branches_name = NULL, min_expr = NULL, 
                                               cell_size = 0.75, norm_method = c("vstExprs", "log"), nrow = NULL, 
                                               ncol = 1, panel_order = NULL, color_by = "Branch", trend_formula = "~sm.ns(Pseudotime, df=3)", 
                                               label_by_short_name = TRUE, TPM = FALSE, cores = 1) {
  if (TPM) {
    exprs(cds) <- esApply(cds, 2, function(x) x/sum(x) * 
                            1e+06)
  }
  if (!(all(branches %in% pData(cds)$State)) & length(branches) == 1) {
    stop("This function only allows to make multiple branch plots where branches is included in the pData")
  }
  branch_label <- branches
  if (!is.null(branches_name)) {
    if (length(branches) != length(branches_name)) {
      stop("branches_name should have the same length as branches")
    }
    branch_label <- branches_name
  }
  g <- cds@minSpanningTree
  m <- NULL
  cds_exprs <- NULL
  for (branch_in in branches) {
    branches_cells <- row.names(subset(pData(cds), State == 
                                         branch_in))
    root_state <- subset(pData(cds), Pseudotime == 0)[, "State"]
    root_state_cells <- row.names(subset(pData(cds), State == 
                                           root_state))
    if (cds@dim_reduce_type != "ICA") {
      root_state_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, 
      ], sep = ""))
      branches_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, 
      ], sep = ""))
    }
    root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 
                                          1)]
    tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 
                                       1)]
    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    if (cds@dim_reduce_type != "ICA") {
      pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
      path_cells <- row.names(pc_ind)[paste("Y_", pc_ind[, 
                                                         1], sep = "") %in% path_cells]
    }
    cds_subset <- cds[, path_cells]
    newdata <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, 
                          row.names = colnames(cds_subset))
    tmp <- t(esApply(cds_subset, 1, function(x) lowess(x[order(pData(cds_subset)$Pseudotime)])$y))
    colnames(tmp) <- colnames(cds_subset)[order(pData(cds_subset)$Pseudotime)]
    cds_exprs_tmp <- reshape2::melt(log2(tmp + 1))
    cds_exprs_tmp <- reshape2::melt(tmp)
    colnames(cds_exprs_tmp) <- c("f_id", "Cell", "expression")
    cds_exprs_tmp$Branch <- branch_label[which(branches == 
                                                 branch_in)]
    if (is.null(cds_exprs)) 
      cds_exprs <- cds_exprs_tmp
    else cds_exprs <- rbind(cds_exprs, cds_exprs_tmp)
    if (is.null(m)) 
      m <- tmp
    else m <- cbind(m, tmp)
  }
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    m = vstExprs(cds, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  if (is.null(min_expr)) {
    min_expr <- cds@lowerDetectionLimit
  }
  cds_pData <- pData(cds)
  cds_fData <- fData(cds)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs <- plyr::ddply(cds_exprs, .(Branch), mutate, Pseudotime = (Pseudotime - 
                                                                         min(Pseudotime)) * 100/(max(Pseudotime) - min(Pseudotime)))
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_line(aes_string(color = color_by), size = I(cell_size))
  }
  q <- q + facet_wrap(~feature_label, nrow = nrow, ncol = ncol, 
                      scales = "free_y")
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  q <- q + monocle_theme_opts()
  q + expand_limits(y = min_expr)
}









