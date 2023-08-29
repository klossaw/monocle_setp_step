library(monocle)
library(monocle3)
library(viridis)
library(tidyverse)

expression_matrix <- readRDS("../expression_matrix2.Rds")
cell_metadata <- readRDS("../cell_metadata2.Rds")
gene_annotation <- readRDS("../gene_annotation2.Rds")
# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)
# cds <- preprocess_cds(cds, num_dim = 50)
# cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
# cds <- reduce_dimension(cds)
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# cds <- cluster_cells(cds)
# plot_cells(cds, color_cells_by = "partition")
# cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "cell.type",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)
# cds <- order_cells(cds)
# mm <- plot_cells(cds,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5)
# # write_rds(cds, file = "monocle3.rds")
# 
# ggplot(mm$data, aes(x = data_dim_1, y = data_dim_2, color = cell_color)) + 
#   geom_point(size = 0.3) + theme_classic() + scale_color_viridis()
# ggplot(mm$data, aes(x = cell_color, y = raw.embryo.time, color = cell_color))  + 
#   geom_point(size = 1,alpha = 0.5) + theme_classic() + scale_color_viridis()
# ggplot(mm$data, aes(x = cell_color, y = embryo.time, color = cell_color)) + 
#   geom_point(size = 1,alpha = 0.5) + theme_classic() + scale_color_viridis()
# ggplot(mm$data, aes(x = cell_color, y = time.point, color = cell_color)) + 
#   geom_point(size = 1,alpha = 0.5) + theme_classic() + scale_color_viridis()
# 
# mm$data$time.point %>% unique()


# cds <- "monocle3.rds" %>% read_rds()
# mm <- plot_cells(cds,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5)
# 
# 
# cell_metadata <- cell_metadata %>% left_join(mm$data)
# rownames(cell_metadata) <- colnames(expression_matrix)
# colnames(cell_metadata)[24] <- "monocle3_time"
# colnames(cell_metadata)[23] <- "monocle3_cluster"
# 
# cell_metadata <- cell_metadata[,-c(20,21,22)]


conflicted::conflicts_prefer(Biobase::fData)
conflicted::conflicts_prefer(Biobase::pData)
conflicted::conflicts_prefer(Biobase::exprs)

# monocle2obj <- monocle::newCellDataSet(expression_matrix,
#                        phenoData = new("AnnotatedDataFrame", data = cell_metadata),
#                        featureData = new("AnnotatedDataFrame", data = gene_annotation),
#                        lowerDetectionLimit = 0.1,
#                        expressionFamily = negbinomial.size())
monocle2obj <- monocle::newCellDataSet(expression_matrix,
                                       phenoData = new("AnnotatedDataFrame", data = cell_metadata[-5]),
                                       featureData = new("AnnotatedDataFrame", data = gene_annotation),
                                       lowerDetectionLimit = 0.1,
                                       expressionFamily = negbinomial.size())
# str(monocle2obj)
exprs(monocle2obj)
fData(monocle2obj)
pData(monocle2obj)

monocle2obj <- estimateSizeFactors(monocle2obj)
monocle2obj <- estimateDispersions(monocle2obj)
monocle2obj <- monocle::detectGenes(monocle2obj, min_expr = 0.05)

sum(fData(monocle2obj)$num_cells_expressed > 300)
sum(pData(monocle2obj)$num_genes_expressed > 200)

ftf <- fData(monocle2obj)$num_cells_expressed > 300
cft <- pData(monocle2obj)$num_genes_expressed > 200
monocle2obj <- monocle2obj[ftf, cft]
monocle2obj <- monocle2obj[,!is.na(monocle2obj$cell.type)]
dim(monocle2obj)

diff_test_resx <- differentialGeneTest(monocle2obj, fullModelFormulaStr = "~ cell.type + bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading", 
                                      reducedModelFormulaStr = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

head(diff_test_resx)

# ========= 差异算法 ==================
cds = monocle2obj
fullModelFormulaStr = "~ cell.type"
reducedModelFormulaStr = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading"
relative_expr = TRUE
cores = 1
verbose = FALSE

status <- NA
all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
pd <- pData(cds)

# diff_test_res <- smartEsApply(cds[1:10,], 1, diff_test_helper,
#                               convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr,
#                               reducedModelFormulaStr = reducedModelFormulaStr,
#                               expressionFamily = cds@expressionFamily, relative_expr = relative_expr,
#                               disp_func = cds@dispFitInfo[["blind"]]$disp_func,
#                               verbose = verbose)
# 
# diff_test_helper <- function (x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, 
#                               relative_expr, weights, disp_func = NULL, verbose = FALSE) {
#   
#   # print(Size_Factor)
#   # print(cell.type)
#   
#   reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, 
#                                   sep = "")
#   fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, 
#                                sep = "")
#   x_orig <- x
#   disp_guess <- 0
#   if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
#     if (relative_expr == TRUE) {
#       x <- x/Size_Factor
#     }
#     f_expression <- round(x)
#     if (is.null(disp_func) == FALSE) {
#       disp_guess <- monocle:::calculate_NB_dispersion_hint(disp_func, 
#                                                  round(x_orig))
#       if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE) {
#         if (expressionFamily@vfamily == "negbinomial") 
#           expressionFamily <- negbinomial(isize = 1/disp_guess)
#         else 
#           expressionFamily <- negbinomial.size(size = 1/disp_guess)
#       }
#     }
#   } else if (expressionFamily@vfamily %in% c("uninormal")) {
#     f_expression <- x
#   } else if (expressionFamily@vfamily %in% c("binomialff")) {
#     f_expression <- x
#   } else {
#     f_expression <- log10(x)
#   }
#   
#   test_res <- tryCatch({
#     if (expressionFamily@vfamily %in% c("binomialff")) {
#       if (verbose) {
#         full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
#                                      epsilon = 0.1, family = expressionFamily)
#         reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
#                                         epsilon = 0.1, family = expressionFamily)
#       } else {
#         
#         full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
#                                                       epsilon = 0.1, family = expressionFamily))
#         reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
#                                                          epsilon = 0.1, family = expressionFamily))
#       }
#     } else {
#       if (verbose) {
#         print(fullModelFormulaStr)
#         full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
#                                      epsilon = 0.1, family = expressionFamily)
#         reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
#                                         epsilon = 0.1, family = expressionFamily)
#       } else {
#         full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
#                                                       epsilon = 0.1, family = expressionFamily))
#         reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
#                                                          epsilon = 0.1, family = expressionFamily))
#       }
#     }
#     
#     compareModels(list(full_model_fit), list(reduced_model_fit))
#   }, error = function(e) {
#     if (verbose) 
#       print(e)
#     data.frame(status = "FAIL", family = expressionFamily@vfamily, 
#                pval = 1, qval = 1)
#   })
#   test_res
# }
head(diff_test_resx,2)
Sp_X <- exprs(cds)[c("WBGene00010957","WBGene00010958"),]
Sp_X <- Matrix::t(Sp_X)
FUN <- diff_test_helper

Size_Factor <- pData(cds)$Size_Factor
cell.type <- pData(cds)$cell.type
bg.300.loading <- pData(cds)$bg.300.loading
bg.400.loading <- pData(cds)$bg.400.loading
bg.500.1.loading <- pData(cds)$bg.500.1.loading
bg.500.2.loading <- pData(cds)$bg.500.2.loading
bg.r17.loading <- pData(cds)$bg.r17.loading
bg.b01.loading <- pData(cds)$bg.b01.loading
bg.b02.loading <- pData(cds)$bg.b02.loading

res <- lapply(colnames(Sp_X), function(i) {
  diff_test_helper(
    as.matrix(Sp_X[, i]), fullModelFormulaStr = fullModelFormulaStr,
      reducedModelFormulaStr = reducedModelFormulaStr,
      expressionFamily = cds@expressionFamily,
      relative_expr = relative_expr,
      disp_func = cds@dispFitInfo[["blind"]]$disp_func,
      verbose = T)
})
names(res) <- colnames(Sp_X)

expressionFamily = cds@expressionFamily
relative_expr = relative_expr
disp_func = cds@dispFitInfo[["blind"]]$disp_func
verbose = verbose
fullModelFormulaStr = "~ cell.type"
reducedModelFormulaStr = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading"
x <- as.matrix(Sp_X[, "WBGene00010957"])
reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep = "")
fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep = "")
x_orig <- x
Size_Factor <- pData(cds)$Size_Factor
disp_guess <- 0
if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
  if (relative_expr == TRUE) {
    x <- x/Size_Factor
  }
  f_expression <- round(x)
  if (is.null(disp_func) == FALSE) {
    disp_guess <- monocle:::calculate_NB_dispersion_hint(disp_func,
                                                         round(x_orig))
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
      full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr),
                                   epsilon = 0.1, family = expressionFamily)
      reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr),
                                      epsilon = 0.1, family = expressionFamily)
    } else {
      full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr),
                                                    epsilon = 0.1, family = expressionFamily))
      reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr),
                                                       epsilon = 0.1, family = expressionFamily))
    }
  } else {
    if (verbose) {
      full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr),
                                   epsilon = 0.1, family = expressionFamily)
      reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr),
                                      epsilon = 0.1, family = expressionFamily)
    } else {
      full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr),
                                                    epsilon = 0.1, family = expressionFamily))
      reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr),
                                                       epsilon = 0.1, family = expressionFamily))
    }
  }
  compareModels(list(full_model_fit), list(reduced_model_fit))
}, error = function(e) {
  if (verbose)
    print(e)
  data.frame(status = "FAIL", family = expressionFamily@vfamily,
             pval = 1, qval = 1)
})
test_res

diff_test_res <- do.call(rbind.data.frame, res)
diff_test_res$qval <- 1
diff_test_res$qval[which(diff_test_res$status == "OK")] <-
  p.adjust(subset(diff_test_res, status == "OK")[, "pval"], method = "BH")
diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
row.names(diff_test_res) <- diff_test_res[, 1]
diff_test_res[, 1] <- NULL
diff_test_res[row.names(cds), ]
# ===========================
# diff_test_res0 <- differentialGeneTest(monocle2obj, fullModelFormulaStr = "~ cell.type")

ordering_genes <- row.names (subset(diff_test_resx, qval < 0.01))
monocle2obj <- setOrderingFilter(monocle2obj, ordering_genes)
plot_ordering_genes(monocle2obj)
monocle2obj <- reduceDimension(monocle2obj, max_components = 2, method = 'DDRTree')

plot(minSpanningTree(monocle2obj), edge.arrow.size=.4,vertex.label=NA, layout = layout_with_fr,vertex.size = 5)

# ========= DDRTree ==================
cds = monocle2obj
max_components = 2
reduction_method = "DDRTree"
norm_method = "log"
pseudo_expr = 1
relative_expr = TRUE
auto_param_selection = TRUE
verbose = TRUE
scaling = TRUE
residualModelFormulaStr = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading"

set.seed(2016)
FM <- monocle:::normalize_expr_data(cds, norm_method, pseudo_expr)
xm <- Matrix::rowMeans(FM)
xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
FM <- FM[xsd > 0, ]


X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),  data = pData(cds), drop.unused.levels = TRUE)
fit <- limma::lmFit(FM, X.model_mat)
beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0
FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])


FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
FM <- FM[!is.na(row.names(FM)), ]
FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ]

ncenter <- monocle:::cal_ncenter(ncol(FM))
cal_ncenter <- function (ncells, ncells_limit = 100) {
  round(2 * ncells_limit * log(ncells)/(log(ncells) + log(ncells_limit)))
}


extra_arguments <- list()
ddr_args <- c(list(X = FM, dimensions = max_components, 
                   ncenter = ncenter, verbose = verbose), 
              extra_arguments[names(extra_arguments) %in% 
                                c("initial_method", "maxIter", "sigma", "lambda", 
                                  "param.gamma", "tol")])
ddrtree_res <- do.call(DDRTree, ddr_args)
  

names(ddr_args)
X = FM
dimensions = ddr_args$dimensions
initial_method = NULL
maxIter = 20
sigma = 0.001
lambda = NULL
ncenter = ddr_args$ncenter
param.gamma = 10
tol = 0.001
verbose = T

D <- nrow(X)
N <- ncol(X)
# 基因的前2个主成分向量
W <- pca_projection_R(X %*% t(X), dimensions)
pca_projection_R <- function (C, L) {
  if (L >= min(dim(C))) {
    eigen_res <- eigen(C)
    U <- eigen_res$vector
    V <- eigen_res$value
    eig_sort <- sort(V, decreasing = T, index.return = T)
    eig_idx <- eig_sort$ix
    W <- U[, eig_idx[1:L]]
    return(W)
  }
  else {
    initial_v <- as.matrix(qnorm(1:(ncol(C) + 1)/(ncol(C) + 
                                                    1))[1:ncol(C)])
    eigen_res <- irlba::irlba(C, nv = L, v = initial_v)
    U <- eigen_res$u
    V <- eigen_res$v
    return(V)
  }
}

# 点迭代初值
Z <- t(W) %*% X

# 抓取骨干细胞
K <- ncenter
centers = t(Z)[seq(1, ncol(Z), length.out = K), ]
kmean_res <- kmeans(t(Z), K, centers = centers)
Y <- kmean_res$centers
Y <- t(Y)

# lambda
lambda = 5 * ncol(X)

ddrtree_res <- DDRTree:::DDRTree_reduce_dim(X, Z, Y, W, dimensions,  maxIter, K, 
                                            sigma, lambda, param.gamma, tol, 
                                            verbose)
names(ddrtree_res)
plot(t(Z))
plot(t(Y))


plot(t(ddrtree_res$Z))
plot(t(ddrtree_res$Y))

colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
colnames(ddrtree_res$Z) <- colnames(FM)

monocle:::reducedDimW(cds) <- ddrtree_res$W
monocle:::reducedDimS(cds) <- ddrtree_res$Z
monocle:::reducedDimK(cds) <- ddrtree_res$Y

cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
adjusted_K <- Matrix::t(reducedDimK(cds))

dp <- as.matrix(dist(adjusted_K))
cellPairwiseDistances(cds) <- dp

library(igraph)
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
minSpanningTree(cds) <- dp_mst
cds@dim_reduce_type <- "DDRTree"
cds <- monocle:::findNearestPointOnMST(cds)
cds


# ===========================
# 第二次运行排序指定根状态
root_cell_candidates <- subset(pData(cds), State == root_state)
dp <- as.matrix(dist(t(reducedDimS(cds)[, row.names(root_cell_candidates)])))
gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp)
tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
diameter <- get.diameter(dp_mst)

root_cell_candidates <- root_cell_candidates[names(diameter),]
if (
  is.null(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) == FALSE && 
  pData(cds)[cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell,]$State == root_state) {
  root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == 
                                                       min(root_cell_candidates$Pseudotime))]
} else {
  root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == 
                                                       max(root_cell_candidates$Pseudotime))]
}
if (length(root_cell) > 1) 
  root_cell <- root_cell[1]
if (cds@dim_reduce_type == "DDRTree") {
  graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell,]
  root_cell = V(minSpanningTree(cds))[graph_point_for_root_cell]$name
}





# 第一次排序不能指定根
conflicted::conflicts_prefer(Biobase::`pData<-`)
cds
root_state = NULL
reverse = NULL


# 第一次排序不能指定根
diameter <- get.diameter(minSpanningTree(cds))
root_cell = names(diameter[1])



dp <- cellPairwiseDistances(cds)
dp_mst <- minSpanningTree(cds)
curr_state <- 1
states = rep(1, ncol(dp))
names(states) <- V(dp_mst)$name

pseudotimes = rep(0, ncol(dp))
names(pseudotimes) <- V(dp_mst)$name

parents = rep(NA, ncol(dp))
names(parents) <- V(dp_mst)$name

mst_traversal <- graph.dfs(dp_mst, root = root_cell, mode = "all", 
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
    curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
    if (degree(dp_mst, v = parent_node_name) > 2) {
      curr_state <- curr_state + 1
    }
  } else {
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
ordering_df


cc_ordering <- ordering_df
pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)),]$pseudo_time
K_old <- reducedDimK(cds)
old_dp <- cellPairwiseDistances(cds)
old_mst <- minSpanningTree(cds)
old_A <- reducedDimA(cds)
old_W <- reducedDimW(cds)



Projection_Method = monocle:::project_point_to_line_segment
dp_mst <- minSpanningTree(cds)
Z <- reducedDimS(cds)
Y <- reducedDimK(cds)
cds <- monocle:::findNearestPointOnMST(cds)
closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
closest_vertex_names <- colnames(Y)[closest_vertex]
closest_vertex_df <- as.matrix(closest_vertex)
row.names(closest_vertex_df) <- row.names(closest_vertex)
tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
for (i in 1:length(closest_vertex)) {
  neighbors <- names(igraph::V(dp_mst)[suppressWarnings(nei(closest_vertex_names[i], mode = "all"))])
  projection <- NULL
  distance <- NULL
  Z_i <- Z[,i]
  for (neighbor in neighbors) {
    if (closest_vertex_names[i] %in% tip_leaves) {
      tmp <- monocle:::projPointOnLine(Z_i, Y[,c(closest_vertex_names[i], neighbor)])
    } else {
      tmp <- Projection_Method(Z_i, Y[,c(closest_vertex_names[i], neighbor)])
    }
    projection <- rbind(projection, tmp)
    distance <- c(distance, dist(rbind(Z_i, tmp)))
  }
  if ("matrix" %in% class(projection)) 
    projection <- as.matrix(projection)
  P[,i] <- projection[which(distance == min(distance))[1],]
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




minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
root_cell_idx <- which(igraph::V(old_mst)$name == root_cell, arr.ind = T)
cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
if (length(cells_mapped_to_graph_root) == 0) {
  cells_mapped_to_graph_root <- root_cell_idx
}
cells_mapped_to_graph_root <-igraph:: V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
tip_leaves <- names(which(igraph::degree(minSpanningTree(cds)) == 1))
root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
if (is.na(root_cell)) {
  root_cell <- select_root_cell(cds, root_state, reverse)
}
cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
cc_ordering_new_pseudotime <- monocle:::extract_ddrtree_ordering(cds, root_cell)
pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)),]$pseudo_time
if (is.null(root_state) == TRUE) {
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  pData(cds)$State <- cc_ordering[closest_vertex[,1], ]$cell_state
}
monocle:::reducedDimK(cds) <- K_old
cellPairwiseDistances(cds) <- old_dp
minSpanningTree(cds) <- old_mst
reducedDimA(cds) <- old_A
monocle:::reducedDimW(cds) <- old_W
mst_branch_nodes <- igraph::V(minSpanningTree(cds))[which(igraph::degree(minSpanningTree(cds)) > 2)]$name
cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
cds


plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "embryo.time")
plot_cell_trajectory(cds, color_by = "time.point")
plot_cell_trajectory(cds, color_by = "monocle3_time")
plot_cell_trajectory(cds, color_by = "monocle3_cluster")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot(cds$Pseudotime,cds$monocle3_time)

# ===========================

monocle2obj <- orderCells(monocle2obj)
plot(t(monocle2obj@reducedDimK))
plot(t(monocle2obj@reducedDimS))

monocle2obj0@reducedDimS
plot(monocle2obj0@reducedDimS, edge.arrow.size=.4,vertex.label=NA, layout = layout_with_fr,vertex.size = 5)

plot_cell_trajectory(monocle2obj, color_by = "Pseudotime")
plot_cell_trajectory(monocle2obj, color_by = "embryo.time")
plot_cell_trajectory(monocle2obj, color_by = "time.point")
plot_cell_trajectory(monocle2obj, color_by = "monocle3_time")
plot_cell_trajectory(monocle2obj, color_by = "monocle3_cluster")
plot_cell_trajectory(monocle2obj, color_by = "State")
plot_cell_trajectory(monocle2obj, color_by = "Pseudotime")
plot(monocle2obj$Pseudotime,monocle2obj$monocle3_time)

ordering_genes <- row.names (subset(diff_test_res0, qval < 0.01))
monocle2obj0 <- setOrderingFilter(monocle2obj, ordering_genes)
plot_ordering_genes(monocle2obj0)
monocle2obj0 <- reduceDimension(monocle2obj0, max_components = 2, method = 'DDRTree')
monocle2obj0 <- orderCells(monocle2obj0)

plot_cell_trajectory(monocle2obj0, color_by = "Pseudotime")
plot_cell_trajectory(monocle2obj0, color_by = "embryo.time")
plot_cell_trajectory(monocle2obj0, color_by = "time.point")
plot_cell_trajectory(monocle2obj0, color_by = "monocle3_time")
plot_cell_trajectory(monocle2obj0, color_by = "monocle3_cluster")
plot_cell_trajectory(monocle2obj0, color_by = "State")
plot_cell_trajectory(monocle2obj0, color_by = "Pseudotime")


plot(monocle2obj0$Pseudotime,monocle2obj0$monocle3_time)

BEAM_res <- BEAM(monocle2obj, branch_point = 1, cores = 1)




plot(monocle2obj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree, layout = 
       layout_as_tree(
         monocle2obj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree, root = c(1, 1),
         rootlevel = c(1, 1)
       ),
     edge.arrow.size=.4,vertex.label=NA,vertex.size = 5)

plot(monocle2obj@minSpanningTree, layout = 
       layout_as_tree(
         monocle2obj@minSpanningTree,
         rootlevel = c(1, 1)
       ),
     edge.arrow.size=.4,vertex.label=NA,vertex.size = 5)

