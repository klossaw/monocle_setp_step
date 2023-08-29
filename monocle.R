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

##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- pbmc@meta.data 
p_data$celltype <- pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
# library(data.table)
# ##读取数据
# data <- fread("fpkm.txt",data.table = F,header = T)
# pd <-  fread("metadata.txt",data.table = F,header = T)
# fd <- fread("gene_annotations.txt",data.table = F,header = T)
# ##创建
# pd <- new("AnnotatedDataFrame", data = pd)
# fd <- new("AnnotatedDataFrame", data = fd)
# HSMM <- newCellDataSet(as.matrix(data),
#                        phenoData = pd, featureData = fd,
#                        expressionFamily = tobit())
# ###如果数据量大，建议转化为稀疏矩阵
# HSMM <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
#                        phenoData = pd, 
#                        featureData = fd,
#                        expressionFamily = tobit())
# =====================================================================
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。


##使用seurat选择的高变基因⚠️
express_genes <- VariableFeatures(pbmc)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
##使用clusters差异表达基因
deg.cluster <- FindAllMarkers(pbmc)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)




#这一步输入的expressed_genes来自于步骤4。
#⚠️⚠️后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(pbmc) 
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
head(cds@phenoData@data)
plot_cell_trajectory(cds,color_by="Pseudotime ", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank())

pdata <- Biobase::pData(cds)
s.cells <- subset(pdata, State=="7") %>% rownames()
# ================================================================
##选择前4个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc

#指定基因
s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
plotc
# =============================================
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
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

diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)

# ===================================
plot_cell_trajectory(cds, color_by = "State")
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2)
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了

#选前100个基因可视化
BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 2, 
                                 num_clusters = 3, show_rownames = T, return_heatmap = T)
#显著差异基因按热图结果排序并保存
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- beam_res[hp.genes, c("gene_short_name", "pval", "qval")]
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "SITX11", "CEBPD", "TYROBP")))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)

