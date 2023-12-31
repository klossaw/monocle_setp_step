# ======================= Clustering and classifying your cells ===========================
library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)


# Load the data
# expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
# cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
# gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))
# 
# saveRDS(expression_matrix, file = 'expression_matrix.Rds')
# saveRDS(cell_metadata, file = 'cell_metadata.Rds')
# saveRDS(gene_annotation, file = 'gene_annotation.Rds')

expression_matrix <- readRDS('expression_matrix.Rds')
cell_metadata <- readRDS('cell_metadata.Rds')
gene_annotation <- readRDS('gene_annotation.Rds')


# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# 各种格式
# cds <- load_cellranger_data("~/Downloads/10x_data")
# cds <- load_mm_data(mat_path = "~/Downloads/matrix.mtx", 
#                     feature_anno_path = "~/Downloads/features.tsv", 
#                     cell_anno_path = "~/Downloads/barcodes.tsv")
# cds <- new_cell_data_set(as(umi_matrix, "sparseMatrix"),
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_metadata)
# 拼接对象
# make a fake second cds object for demonstration
# cds2 <- cds[1:100,]
# 
# big_cds <- combine_cds(list(cds, cds2))

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")


plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)
cds = align_cds(cds, num_dim = 100, alignment_group = "plate")
cds = reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)

marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)


colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="Germline",
                                                "2"="Body wall muscle",
                                                "3"="Unclassified neurons",
                                                "4"="Vulval precursors",
                                                "5"="Failed QC",
                                                "6"="Seam cells",
                                                "7"="Pharyngeal epithelia",
                                                "8"="Coelomocytes",
                                                "9"="Am/PH sheath cells",
                                                "10"="Failed QC",
                                                "11"="Touch receptor neurons",
                                                "12"="Intestinal/rectal muscle",
                                                "13"="Pharyngeal neurons",
                                                "14"="NA",
                                                "15"="flp-1(+) interneurons",
                                                "16"="Canal associated neurons",
                                                "17"="Ciliated sensory neurons",
                                                "18"="Other interneurons",
                                                "19"="Pharyngeal gland",
                                                "20"="Failed QC",
                                                "21"="Ciliated sensory neurons",
                                                "22"="Oxygen sensory neurons",
                                                "23"="Ciliated sensory neurons",
                                                "24"="Ciliated sensory neurons",
                                                "25"="Ciliated sensory neurons",
                                                "26"="Ciliated sensory neurons",
                                                "27"="Oxygen sensory neurons",
                                                "28"="Ciliated sensory neurons",
                                                "29"="Unclassified neurons",
                                                "30"="Socket cells",
                                                "31"="Failed QC",
                                                "32"="Pharyngeal gland",
                                                "33"="Ciliated sensory neurons",
                                                "34"="Ciliated sensory neurons",
                                                "35"="Ciliated sensory neurons",
                                                "36"="Failed QC",
                                                "37"="Ciliated sensory neurons",
                                                "38"="Pharyngeal muscle")
plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")

# 取子集
cds_subset <- choose_cells(cds)
pr_graph_test_res <- graph_test(cds_subset, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-3)
plot_cells(cds_subset, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
cds_subset = cluster_cells(cds_subset, resolution=1e-2)
plot_cells(cds_subset, color_cells_by="cluster")


colData(cds_subset)$assigned_cell_type <- as.character(clusters(cds_subset)[colnames(cds_subset)])
colData(cds_subset)$assigned_cell_type <- dplyr::recode(colData(cds_subset)$assigned_cell_type,
                                                        "1"="Somatic gonad precursors",
                                                        "2"="Somatic gonad precursors",
                                                        "3"="Vulval precursors",
                                                        "4"="Sex myoblasts",
                                                        "5"="Sex myoblasts",
                                                        "6"="Vulval precursors",
                                                        "7"="Failed QC",
                                                        "8"="Vulval precursors",
                                                        "10"="Unclassified neurons",
                                                        "11"="Distal tip cells")
plot_cells(cds_subset, group_cells_by="cluster", color_cells_by="assigned_cell_type")
colData(cds)[colnames(cds_subset),]$assigned_cell_type <- colData(cds_subset)$assigned_cell_type
cds <- cds[,colData(cds)$assigned_cell_type != "Failed QC" | is.na(colData(cds)$assigned_cell_type )]
plot_cells(cds, group_cells_by="partition", 
           color_cells_by="assigned_cell_type", 
           labels_per_group=5)




assigned_type_marker_test_res <- top_markers(cds,
                                             group_cells_by="assigned_cell_type",
                                             reference_cells=1000,
                                             cores=5)

# Require that markers have at least JS specificty score > 0.5 and
# be significant in the logistic test for identifying their cell type:
garnett_markers <- assigned_type_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)
# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% 
  group_by(gene_short_name) %>%
  filter(n() == 1)
generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")



## Install the monocle3 branch of garnett
# BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
# devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
library(garnett)
# install gene database for worm
# BiocManager::install("org.Ce.eg.db")

colData(cds)$garnett_cluster <- clusters(cds)
worm_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = "./marker_file.txt", 
                                         db=org.Ce.eg.db::org.Ce.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         cores=8)
cds <- classify_cells(cds, worm_classifier,
                      db = org.Ce.eg.db::org.Ce.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "ENSEMBL")

plot_cells(cds,
           group_cells_by="partition",
           color_cells_by="cluster_ext_type")

# ceWhole <- readRDS(url("https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole_20191017.RDS"))
# saveRDS(ceWhole, file = 'ceWhole.Rds')
ceWhole <- readRDS('ceWhole.Rds')
cds <- classify_cells(cds, ceWhole,
                      db = org.Ce.eg.db::org.Ce.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "ENSEMBL")

# =================== Constructing single-cell trajectories ===============================
library(monocle3)
library(ggplot2)
library(dplyr)
# expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
# cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
# gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))
# 
# saveRDS(expression_matrix,file = 'expression_matrix2.Rds')
# saveRDS(cell_metadata,file = 'cell_metadata2.Rds')
# saveRDS(gene_annotation,file = 'gene_annotation2.Rds')

expression_matrix <- readRDS('expression_matrix2.Rds')
cell_metadata <- readRDS('cell_metadata2.Rds')
gene_annotation <- readRDS('gene_annotation2.Rds')

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", 
                 residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE, 
           color_cells_by = "cell.type")

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# 3D
cds_sub <- choose_graph_segments(cds)

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
cds_3d_plot_obj

# ================== Differential expression analysis ===========================
ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")
cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")
emb_time_terms %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs %>% filter(term != "(Intercept)") %>%
  dplyr::select(gene_short_name, term, q_value, estimate)

evaluate_fits(gene_fits)
time_batch_models <- fit_models(cds_subset,
                                model_formula_str = "~embryo.time + batch",
                                expression_family="negbinomial")
time_models <- fit_models(cds_subset,
                          model_formula_str = "~embryo.time",
                          expression_family="negbinomial")
View(time_models$model_summary$WBGene00001820)
compare_models(time_batch_models, time_models) %>% 
  dplyr::select(gene_short_name, q_value)



# reload and reprocess the data as described in the 'Clustering and classifying your cells' section
expression_matrix <- readRDS('expression_matrix.Rds')
cell_metadata <- readRDS('cell_metadata.Rds')
gene_annotation <- readRDS('gene_annotation.Rds')
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-5)
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Germline",
                                                 "2"="Body wall muscle",
                                                 "3"="Unclassified neurons",
                                                 "4"="Vulval precursors",
                                                 "5"="Failed QC",
                                                 "6"="Seam cells",
                                                 "7"="Pharyngeal epithelia",
                                                 "8"="Coelomocytes",
                                                 "9"="Am/PH sheath cells",
                                                 "10"="Failed QC",
                                                 "11"="Touch receptor neurons",
                                                 "12"="Intestinal/rectal muscle",
                                                 "13"="Pharyngeal neurons",
                                                 "14"="NA",
                                                 "15"="flp-1(+) interneurons",
                                                 "16"="Canal associated neurons",
                                                 "17"="Ciliated sensory neurons",
                                                 "18"="Other interneurons",
                                                 "19"="Pharyngeal gland",
                                                 "20"="Failed QC",
                                                 "21"="Ciliated sensory neurons",
                                                 "22"="Oxygen sensory neurons",
                                                 "23"="Ciliated sensory neurons",
                                                 "24"="Ciliated sensory neurons",
                                                 "25"="Ciliated sensory neurons",
                                                 "26"="Ciliated sensory neurons",
                                                 "27"="Oxygen sensory neurons",
                                                 "28"="Ciliated sensory neurons",
                                                 "29"="Unclassified neurons",
                                                 "30"="Socket cells",
                                                 "31"="Failed QC",
                                                 "32"="Pharyngeal gland",
                                                 "33"="Ciliated sensory neurons",
                                                 "34"="Ciliated sensory neurons",
                                                 "35"="Ciliated sensory neurons",
                                                 "36"="Failed QC",
                                                 "37"="Ciliated sensory neurons",
                                                 "38"="Pharyngeal muscle")
neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
plot_cells(neurons_cds, color_cells_by="partition")


pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)), 
                                cell_group=partitions(cds)[colnames(neurons_cds)])
agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
plot_cells(neurons_cds, 
           genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)




expression_matrix <- readRDS('expression_matrix2.Rds')
cell_metadata <- readRDS('cell_metadata2.Rds')
gene_annotation <- readRDS('gene_annotation2.Rds')
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds)

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)

cds_subset <- choose_cells(cds)
subset_pr_test_res <- graph_test(cds_subset, 
                                 neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=T,
           show_trajectory_graph=T)









