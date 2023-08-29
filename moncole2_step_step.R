# ========================= input ===============================
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork) #用来做拼图的包，后面的p1|p1|p3在一张图上展示三个图就是这个包的功劳

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc[['cell_type']] <- pbmc@active.ident #将注释结果添加到metadata
levels(pbmc)

# ================================================================
library(monocle)
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
p_data$celltype <- pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",
                             cores=1) 
head(diff)

deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5)
plot_cell_trajectory(cds,color_by="Pseudotime ", size=1,show_backbone=TRUE) 
plot_complex_cell_trajectory(cds,
                             color_by = "State",
                             root_states = 3)


pdata <- Biobase::pData(cds)
s.cells <- subset(pdata, State=="7") %>% rownames()
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")

#指定基因
s.genes <- c("SELL","CCR7","IL7R", 
             "CD84","CCL5","S100A4")
plot_genes_jitter(cds[s.genes[1],], grouping = "State", 
                  color_by = "State")
plot_genes_violin(cds[s.genes,], grouping = "State", 
                  color_by = "State")


Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>%
  pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, 
                            show_rownames=T, return_heatmap=T)
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))
plot_cell_trajectory(cds, color_by = "State")
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了

# ============================= plot_pseudotime_heatmap ========================
plot_pseudotime_heatmap <- function (cds_subset, cluster_rows = TRUE, 
                                     hclust_method = "ward.D2", 
          num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
          add_annotation_col = NULL, show_rownames = FALSE, 
          use_gene_short_name = TRUE, 
          norm_method = c("log", "vstExprs"), scale_max = 3, 
          scale_min = -3, trend_formula = "~sm.ns(Pseudotime, df=3)", 
          return_heatmap = FALSE, cores = 1) 
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  # m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
  #                      relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    # m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    # m = log10(m + pseudocount)
  }
  # m = m[!apply(m, 1, sd) == 0, ]
  # m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  # m = m[is.na(row.names(m)) == FALSE, ]
  # m[is.nan(m)] = 0
  # m[m > scale_max] = scale_max
  # m[m < scale_min] = scale_min
  # heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  # ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
  #                cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
  #                clustering_distance_rows = row_dist, clustering_method = hclust_method, 
  #                cutree_rows = num_clusters, silent = TRUE, filename = NA, 
  #                breaks = bks, border_color = NA, color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                         num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
    ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if (!is.null(annotation_row)) 
      row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row)) 
    row.names(annotation_row) <- row_ann_labels
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  # ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
  #                    cluster_rows = cluster_rows, show_rownames = show_rownames, 
  #                    show_colnames = F, clustering_distance_rows = row_dist, 
  #                    clustering_method = hclust_method, cutree_rows = num_clusters, 
  #                    annotation_row = annotation_row, annotation_col = annotation_col, 
  #                    treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
  #                    border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(ph_res)
  }
}
# =======~~~~~~~~~~   自定义  ~~~~~~~~~~ =======
trend_formula = "~sm.ns(Pseudotime, df=3)"
cds_subset <- cds
pseudocount <- 1
newdata <- data.frame(Pseudotime = seq(
  min(pData(cds_subset)$Pseudotime), 
  max(pData(cds_subset)$Pseudotime), length.out = 100))
m <- genSmoothCurves(cds_subset, cores = 4, 
                     trend_formula = trend_formula,
                     relative_expr = T, new_data = newdata)
m = m[!apply(m, 1, sum) == 0, ]
m = vstExprs(cds_subset, expr_matrix = m)
m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0

scale_max = 3
scale_min = -3
m[m > scale_max] = scale_max
m[m < scale_min] = scale_min
heatmap_matrix <- m

name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
Heatmap(heatmap_matrix, cluster_columns = F, 
        show_row_names = F,
        show_column_names = F,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2), 
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+ 
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name), 
                                 labels = name, 
                                 labels_gp = gpar(fontsize = 10)))
  


# ============================= plot_genes_branched_heatmap ====================
plot_genes_branched_heatmap <- function (cds_subset, branch_point = 1, branch_states = NULL, 
          branch_labels = c("Cell fate 1", "Cell fate 2"), 
          cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 6, 
          hmcols = NULL, branch_colors = c("#979797", "#F05662", 
                                           "#7990C8"), add_annotation_row = NULL, add_annotation_col = NULL, 
          show_rownames = FALSE, use_gene_short_name = TRUE, scale_max = 3, 
          scale_min = -3, norm_method = c("log", "vstExprs"), 
          trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", 
          return_heatmap = FALSE, cores = 1, ...) 
{
  cds <- NA
  # new_cds <- buildBranchCellDataSet(cds_subset, branch_states = branch_states, 
  #                                   branch_point = branch_point, progenitor_method = "duplicate", 
  #                                   ...)
  # new_cds@dispFitInfo <- cds_subset@dispFitInfo
  # if (is.null(branch_states)) {
  #   progenitor_state <- subset(pData(cds_subset), Pseudotime == 
  #                                0)[, "State"]
  #   branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  # }
  # col_gap_ind <- 101
  # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
  #                        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
  # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100), 
  #                        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  # BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores = cores, 
  #                                   trend_formula = trend_formula, relative_expr = T, new_data = rbind(newdataA, 
  #                                                                                                      newdataB))
  # BranchA_exprs <- BranchAB_exprs[, 1:100]
  # BranchB_exprs <- BranchAB_exprs[, 101:200]
  # common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == 
  #                                                     setdiff(pData(new_cds)$State, 
  #                                                             branch_states), ])
  # BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 
  #                                                "Pseudotime"])))
  # BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 
  #                                         "Pseudotime"]))
  # BranchB_num <- BranchA_num
  # norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs") {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
  }
  else if (norm_method == "log") {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  # heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], 
  #                         BranchB_exprs)
  # heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, 
  #                                        sd) == 0, ]
  # heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), 
  #                                  center = TRUE))
  # heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == 
  #                                   FALSE, ]
  # heatmap_matrix[is.nan(heatmap_matrix)] = 0
  # heatmap_matrix[heatmap_matrix > scale_max] = scale_max
  # heatmap_matrix[heatmap_matrix < scale_min] = scale_min
  # heatmap_matrix_ori <- heatmap_matrix
  # heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 
  #                                                           1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
  # row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  # row_dist[is.na(row_dist)] <- 1
  # exp_rng <- range(heatmap_matrix)
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
  if (is.null(hmcols)) {
    hmcols <- blue2green2red(length(bks) - 1)
  }
  # ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
  #                cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
  #                clustering_distance_rows = row_dist, clustering_method = hclust_method, 
  #                cutree_rows = num_clusters, silent = TRUE, filename = NA, 
  #                breaks = bks, color = hmcols)
  annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                       num_clusters)))
  if (!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
    ])
  }
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), 
                               `Cell Type` = c(rep(branch_labels[1], BranchA_num), 
                                               rep("Pre-branch", 2 * BranchP_num), rep(branch_labels[2], 
                                                                                       BranchB_num)))
  colnames(annotation_col) <- "Cell Type"
  if (!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), 
    ])$gene_short_name, 1])
  }
  names(branch_colors) <- c("Pre-branch", branch_labels[1], 
                            branch_labels[2])
  annotation_colors = list(`Cell Type` = branch_colors)
  names(annotation_colors$`Cell Type`) = c("Pre-branch", 
                                           branch_labels)
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  # ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
  #                    cluster_rows = TRUE, show_rownames = show_rownames, show_colnames = F, 
  #                    clustering_distance_rows = row_dist, clustering_method = hclust_method, 
  #                    cutree_rows = num_clusters, annotation_row = annotation_row, 
  #                    annotation_col = annotation_col, annotation_colors = annotation_colors, 
  #                    gaps_col = col_gap_ind, treeheight_row = 20, breaks = bks, 
  #                    fontsize = 6, color = hmcols, border_color = NA, silent = TRUE)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, 
                heatmap_matrix = heatmap_matrix, heatmap_matrix_ori = heatmap_matrix_ori, 
                ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, 
                hmcols = hmcols, annotation_colors = annotation_colors, 
                annotation_row = annotation_row, annotation_col = annotation_col, 
                ph_res = ph_res))
  }
}
# =======~~~~~~~~~~~~~~~~~~~~ 自定义  ~~~~~~~~~~~~~~~~~~~~~~~~=======
branch_labels = c("Cell fate 1", "Cell fate 2")  
cds_subset <- cds[row.names(subset(BEAM_res,
                                   qval < 1e-4)),]
new_cds <- buildBranchCellDataSet(cds_subset, 
                                  branch_point = 1, 
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula, 
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
    
    
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], 
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, 
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), 
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) == 
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), 
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num), 
                                             rep("Pre-branch", 2 * BranchP_num), 
                                             rep(branch_labels[2], 
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100), 
            rep('Cell fate 2', times = 100))
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[1], 
                          "Cell fate 2" = branch_colors[2],
                          "Pre-branch" = branch_colors[3])),
  na_col = "grey")
Heatmap(heatmap_matrix, cluster_columns = F, 
        show_row_names = F,
        show_column_names = F,
        column_split = csplit,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2), 
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+ 
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name), 
                                 labels = name, 
                                 labels_gp = gpar(fontsize = 10)))
  
# ============================= plot_genes_jitter ==============================

# ============================= plot_genes_violin ==============================

# ============================= plot_genes_in_pseudotime =======================

# ============================= plot_cell_trajectory ===========================

# ============================= plot_complex_cell_trajectory ===================
function (cds, x = 1, y = 2, root_states = NULL, color_by = "State", 
          show_tree = TRUE, show_backbone = TRUE, backbone_color = "black", 
          markers = NULL, show_cell_names = FALSE, cell_size = 1.5, 
          cell_link_size = 0.75, cell_name_size = 2, show_branch_points = TRUE, 
          ...) 
{
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
  }
  else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", 
                                      "SGL-tree")) {
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  if (is.null(reduced_dim_coords)) {
    stop("You must first call reduceDimension() before using this function")
  }
  dp_mst <- minSpanningTree(cds)
  if (is.null(root_states)) {
    if (is.null(lib_info_with_pseudo$Pseudotime)) {
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 
                                                     1][1]
    }
    else root_cell <- row.names(subset(lib_info_with_pseudo, 
                                       Pseudotime == 0))
    if (cds@dim_reduce_type != "ICA") 
      root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, 
      ]]
  }
  else {
    candidate_root_cells <- row.names(subset(pData(cds), 
                                             State %in% root_states))
    if (cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, 
                                                     candidate_root_cells) == 1)]
    }
    else {
      Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, 
      ]]
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, 
                                                       Y_candidate_root_cells) == 1)]
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
  edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                             "prin_graph_dim_1", "prin_graph_dim_2")], 
                   by.x = "target", by.y = "sample_name", all = TRUE)
  edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "target_prin_graph_dim_2"))
  if (cds@dim_reduce_type == "ICA") {
    S_matrix <- tree_coords[, ]
  }
  else if (cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", 
                                      "SGL-tree")) {
    S_matrix <- tree_coords[closest_vertex, ]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "row.names")
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), 
      ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", 
                                        "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, 
                                    I(cell_size))) + facet_wrap(~feature_label)
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    if (class(data_df[, color_by]) == "numeric") {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", 
                                                     color_by, " + 0.1)")), size = I(cell_size), 
                           na.rm = TRUE, height = 5) + scale_color_viridis(name = paste0("log10(", 
                                                                                         color_by, ")"), ...)
    }
    else {
      g <- g + geom_jitter(aes_string(color = color_by), 
                           size = I(cell_size), na.rm = TRUE, height = 5)
    }
  }
  else {
    if (class(data_df[, color_by]) == "numeric") {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", 
                                                     color_by, " + 0.1)")), size = I(cell_size), 
                           na.rm = TRUE, height = 5) + scale_color_viridis(name = paste0("log10(", 
                                                                                         color_by, " + 0.1)"), ...)
    }
    else {
      g <- g + geom_jitter(aes_string(color = color_by), 
                           size = I(cell_size), na.rm = TRUE, height = 5)
    }
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[, 
                                                                          c("sample_name", "source_prin_graph_dim_1", 
                                                                            "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, 
                                              mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), 
    ]
    g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2"), size = 2 * 
                          cell_size, na.rm = TRUE, data = branch_point_df) + 
      geom_text(aes_string(x = "source_prin_graph_dim_1", 
                           y = "source_prin_graph_dim_2", label = "branch_point_idx"), 
                size = 1.5 * cell_size, color = "white", 
                na.rm = TRUE, data = branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  g <- g + theme(strip.background = element_rect(colour = "white", 
                                                 fill = "white")) + theme(panel.border = element_blank()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(legend.key = element_blank()) + xlab("") + 
    ylab("") + theme(legend.position = "top", 
                     legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(line = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  g
}
# =================== ~~~~~~~~~ 自定义 ~~~~~~~~~~~ ==============================
library(igraph)
cds=cds
x = 1
y = 2
root_states = NULL
color_by = "State"
show_tree = TRUE
show_backbone = TRUE
backbone_color = "black"
markers = NULL
show_cell_names = FALSE
cell_size = 1.5
cell_link_size = 0.75
cell_name_size = 2
show_branch_points = TRUE
root_states = 3

gene_short_name <- NA
sample_name <- NA
data_dim_1 <- NA
data_dim_2 <- NA
lib_info_with_pseudo <- pData(cds)
  
if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", 
                                      "SGL-tree")) {
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
}

dp_mst <- minSpanningTree(cds)
# root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))
# candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
# 默认方法
Y_candidate_root_cells <- V(dp_mst)$name[
  cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[
    candidate_root_cells,]]
root_cell <- Y_candidate_root_cells[which(degree(dp_mst,
                                                 Y_candidate_root_cells) == 1)]

tree_coords <- layout_as_tree(dp_mst, root = root_cell)
ica_space_df <- data.frame(tree_coords)
row.names(ica_space_df) <- colnames(reduced_dim_coords)
colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c("source", "target")
edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                   by.y = "source", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "source_prin_graph_dim_2"))
edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                             "prin_graph_dim_1", "prin_graph_dim_2")], 
                   by.x = "target", by.y = "sample_name", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                     prin_graph_dim_2 = "target_prin_graph_dim_2"))
S_matrix <- tree_coords[closest_vertex, ]
closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  
data_df <- data.frame(S_matrix)
row.names(data_df) <- colnames(reducedDimS(cds))
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- row.names(data_df)
data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name",
                 by.y = "row.names")
  
markers_exprs <- NULL


g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))


g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = cell_link_size, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  


g <- g + geom_jitter(aes_string(color = color_by), 
                     size = I(cell_size), na.rm = TRUE, height = 5)
    
  
mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[, 
                                                                          c("sample_name", "source_prin_graph_dim_1", 
                                                                            "source_prin_graph_dim_2")]
branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, 
                                              mst_branch_nodes)
branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), 
    ]
g <- g + geom_point(aes_string(x = "source_prin_graph_dim_1", 
                                   y = "source_prin_graph_dim_2"), size = 2 * 
                          cell_size, na.rm = TRUE, data = branch_point_df) + 
      geom_text(aes_string(x = "source_prin_graph_dim_1", 
                           y = "source_prin_graph_dim_2", label = "branch_point_idx"), 
                size = 1.5 * cell_size, color = "white", 
                na.rm = TRUE, data = branch_point_df)

g <- g + theme(strip.background = element_rect(colour = "white", 
                                                 fill = "white")) + theme(panel.border = element_blank()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(legend.key = element_blank()) + xlab("") + 
    ylab("") + theme(legend.position = "top", 
                     legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(line = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  g

# ============================= setOrderingFilter ==============================
fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
# ============================= reduceDimension ================================
  function (cds, max_components = 2, reduction_method = c("DDRTree", 
                                                          "ICA", "tSNE", "SimplePPT", "L1-graph", 
                                                          "SGL-tree"), norm_method = c("log", "vstExprs", 
                                                                                       "none"), residualModelFormulaStr = NULL, pseudo_expr = 1, 
            relative_expr = TRUE, auto_param_selection = TRUE, verbose = FALSE, 
            scaling = TRUE, ...) 
  {
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
    }
    else {
      X.model_mat <- NULL
    }
    if (scaling) {
      FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      FM <- FM[!is.na(row.names(FM)), ]
    }
    else FM <- as.matrix(FM)
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
      gp <- graph.adjacency(dp, mode = "undirected", 
                            weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "function_passed"
    }
    else {
      reduction_method <- match.arg(reduction_method)
      if (reduction_method == "tSNE") {
        if (verbose) 
          message("Remove noise by PCA ...")
        if ("num_dim" %in% names(extra_arguments)) {
          num_dim <- extra_arguments$num_dim
        }
        else {
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
      }
      else if (reduction_method == "ICA") {
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
        gp <- graph.adjacency(dp, mode = "undirected", 
                              weighted = TRUE)
        dp_mst <- minimum.spanning.tree(gp)
        minSpanningTree(cds) <- dp_mst
        cds@dim_reduce_type <- "ICA"
      }
      else if (reduction_method == "DDRTree") {
        if (verbose) 
          message("Learning principal graph with DDRTree")
        if (auto_param_selection & ncol(cds) >= 100) {
          if ("ncenter" %in% names(extra_arguments)) 
            ncenter <- extra_arguments$ncenter
          else ncenter <- cal_ncenter(ncol(FM))
          ddr_args <- c(list(X = FM, dimensions = max_components, 
                             ncenter = ncenter, verbose = verbose), extra_arguments[names(extra_arguments) %in% 
                                                                                      c("initial_method", "maxIter", 
                                                                                        "sigma", "lambda", "param.gamma", 
                                                                                        "tol")])
          ddrtree_res <- do.call(DDRTree, ddr_args)
        }
        else {
          ddrtree_res <- DDRTree(FM, max_components, verbose = verbose, 
                                 ...)
        }
        if (ncol(ddrtree_res$Y) == ncol(cds)) 
          colnames(ddrtree_res$Y) <- colnames(FM)
        else colnames(ddrtree_res$Y) <- paste("Y_", 
                                              1:ncol(ddrtree_res$Y), sep = "")
        colnames(ddrtree_res$Z) <- colnames(FM)
        reducedDimW(cds) <- ddrtree_res$W
        reducedDimS(cds) <- ddrtree_res$Z
        reducedDimK(cds) <- ddrtree_res$Y
        cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals
        adjusted_K <- Matrix::t(reducedDimK(cds))
        dp <- as.matrix(dist(adjusted_K))
        cellPairwiseDistances(cds) <- dp
        gp <- graph.adjacency(dp, mode = "undirected", 
                              weighted = TRUE)
        dp_mst <- minimum.spanning.tree(gp)
        minSpanningTree(cds) <- dp_mst
        cds@dim_reduce_type <- "DDRTree"
        cds <- findNearestPointOnMST(cds)
      }
      else {
        stop("Error: unrecognized dimensionality reduction method")
      }
    }
    cds
  }
# ============================= orderCells =====================================
  function (cds, root_state = NULL, num_paths = NULL, reverse = NULL) 
  {
    if (class(cds)[1] != "CellDataSet") {
      stop("Error cds is not of type 'CellDataSet'")
    }
    if (is.null(cds@dim_reduce_type)) {
      stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }
    if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 
              0))) {
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
      gp <- graph.adjacency(dp, mode = "undirected", 
                            weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- dp_mst
      next_node <<- 0
      res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
      cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
      order_list <- extract_good_branched_ordering(res$subtree, 
                                                   res$root, cellPairwiseDistances(cds), num_paths, 
                                                   FALSE)
      cc_ordering <- order_list$ordering_df
      row.names(cc_ordering) <- cc_ordering$sample_name
      minSpanningTree(cds) <- as.undirected(order_list$cell_ordering_tree)
      pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
      ]$pseudo_time
      pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
      ]$cell_state
      mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                          2)]$name
      minSpanningTree(cds) <- dp_mst
      cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- as.undirected(order_list$cell_ordering_tree)
    }
    else if (cds@dim_reduce_type == "DDRTree") {
      if (is.null(num_paths) == FALSE) {
        message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
      }
      cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
      pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
      ]$pseudo_time
      K_old <- reducedDimK(cds)
      old_dp <- cellPairwiseDistances(cds)
      old_mst <- minSpanningTree(cds)
      old_A <- reducedDimA(cds)
      old_W <- reducedDimW(cds)
      cds <- project2MST(cds, project_point_to_line_segment)
      minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
      root_cell_idx <- which(V(old_mst)$name == root_cell, 
                             arr.ind = T)
      cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == 
                                            root_cell_idx)
      if (length(cells_mapped_to_graph_root) == 0) {
        cells_mapped_to_graph_root <- root_cell_idx
      }
      cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
      tip_leaves <- names(which(degree(minSpanningTree(cds)) == 
                                  1))
      root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% 
                                                tip_leaves][1]
      if (is.na(root_cell)) {
        root_cell <- select_root_cell(cds, root_state, reverse)
      }
      cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
      cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, 
                                                             root_cell)
      pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), 
      ]$pseudo_time
      if (is.null(root_state) == TRUE) {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        pData(cds)$State <- cc_ordering[closest_vertex[, 
                                                       1], ]$cell_state
      }
      reducedDimK(cds) <- K_old
      cellPairwiseDistances(cds) <- old_dp
      minSpanningTree(cds) <- old_mst
      reducedDimA(cds) <- old_A
      reducedDimW(cds) <- old_W
      mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                          2)]$name
    }
    else if (cds@dim_reduce_type == "SimplePPT") {
      if (is.null(num_paths) == FALSE) {
        message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
      }
      cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
      pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
      ]$pseudo_time
      pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
      ]$cell_state
      mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                          2)]$name
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
    cds
  }
# ============================= BEAM ===========================================

# ============================= differentialGeneTest ===========================================
  function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
            reducedModelFormulaStr = "~1", relative_expr = TRUE, 
            cores = 1, verbose = FALSE) 
  {
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
    if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", 
                                                             "negbinomial.size")) {
      if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
        stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
      }
    }
    if (cores > 1) {
      diff_test_res <- mcesApply(cds, 1, diff_test_helper, 
                                 c("BiocGenerics", "VGAM", "Matrix"), 
                                 cores = cores, fullModelFormulaStr = fullModelFormulaStr, 
                                 reducedModelFormulaStr = reducedModelFormulaStr, 
                                 expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                                 disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                                 verbose = verbose)
      diff_test_res
    }
    else {
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
    diff_test_res$qval[which(diff_test_res$status == "OK")] <- p.adjust(subset(diff_test_res, 
                                                                               status == "OK")[, "pval"], method = "BH")
    diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
    row.names(diff_test_res) <- diff_test_res[, 1]
    diff_test_res[, 1] <- NULL
    diff_test_res[row.names(cds), ]
  }
# ============================= genSmoothCurves ===========================================
# ============================= vstExprs ===========================================

######################


