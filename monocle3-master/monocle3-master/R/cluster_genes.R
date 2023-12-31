#' Cluster genes into modules that are co-expressed across cells.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param reduction_method The dimensionality reduction method used to generate the lower dimensional space in which genes will be clustered. Currently only UMAP is supported.
#' @param max_components The number of dimensions in which to cluster genes into modules.
#' @param umap.metric Metric used by UMAP for measuring similarity between genes .
#' @param umap.min_dist Minimum distance parameter passed to UMAP.
#' @param umap.n_neighbors Number of nearest neighbors used by UMAP.
#' @param umap.fast_sgd Whether to allow UMAP to perform fast stochastic gradient descent. Defaults to TRUE. Setting FALSE will result in slower, but deterministic behavior (if cores=1).
#' @param umap.nn_method The method used for nearest neighbor network construction during UMAP.
#' @param k number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa. Default to be 20
#' @param leiden_iter Integer number of iterations used for Leiden clustering. The clustering result with the largest modularity score is used as the final clustering result.  Default to be 1.
#' @param partition_qval Significance threshold used in Louvain community graph partitioning.
#' @param weight A logic argument to determine whether or not we will use
#'   Jaccard coefficient for two nearest neighbors (based on the overlapping of
#'   their kNN) as the weight used for Louvain clustering. Default to be FALSE.
#' @param resolution Resolution parameter passed to Louvain. Can be a list. If
#'   so, this method will evaluate modularity at each resolution and use the
#'   one with the highest value.
#' @param random_seed  the seed used by the random number generator in Leiden.
#' @param cores number of cores computer should use to execute function
#' @param verbose Whether or not verbose output is printed.
#' @param preprocess_method a string specifying the low-dimensional space
#'   to use for gene loadings, currently either PCA or LSI. Default is
#'   "PCA".
#' @param nn_control An optional list of parameters used to make the nearest
#'  neighbor index. See the set_nn_control help for detailed information.
#' @param ... Additional arguments passed to UMAP and Louvain analysis.
#'
#' @return A dataframe with genes and the modules to which they are assigned.
#'
#' @examples
#'   \dontrun{
#'      expression_matrix <- readRDS(system.file('extdata',
#'                                               'worm_l2/worm_l2_expression_matrix.rds',
#'                                               package='monocle3'))
#'      cell_metadata <- readRDS(system.file('extdata',
#'                                           'worm_l2/worm_l2_coldata.rds',
#'                                           package='monocle3'))
#'      gene_metadata <- readRDS(system.file('extdata',
#'                                           'worm_l2/worm_l2_rowdata.rds',
#'                                           package='monocle3'))
#'
#'      cds <- new_cell_data_set(expression_data=expression_matrix,
#'                               cell_metadata=cell_metadata,
#'                               gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds, num_dim = 100)
#'     cds <- reduce_dimension(cds)
#'     cds <- cluster_cells(cds, resolution=1e-5)
#'     colData(cds)$assigned_cell_type <- as.character(partitions(cds))
#'     colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
#'                                                     "1"="Germline",
#'                                                     "2"="Body wall muscle",
#'                                                     "3"="Unclassified neurons",
#'                                                     "4"="Vulval precursors",
#'                                                     "5"="Failed QC",
#'                                                     "6"="Seam cells",
#'                                                     "7"="Pharyngeal epithelia",
#'                                                     "8"="Coelomocytes",
#'                                                     "9"="Am/PH sheath cells",
#'                                                     "10"="Failed QC",
#'                                                     "11"="Touch receptor neurons",
#'                                                     "12"="Intestinal/rectal muscle",
#'                                                     "13"="Pharyngeal neurons",
#'                                                     "14"="NA",
#'                                                     "15"="flp-1(+) interneurons",
#'                                                     "16"="Canal associated neurons",
#'                                                     "17"="Ciliated sensory neurons",
#'                                                     "18"="Other interneurons",
#'                                                     "19"="Pharyngeal gland",
#'                                                     "20"="Failed QC",
#'                                                     "21"="Ciliated sensory neurons",
#'                                                     "22"="Oxygen sensory neurons",
#'                                                     "23"="Ciliated sensory neurons",
#'                                                     "24"="Ciliated sensory neurons",
#'                                                     "25"="Ciliated sensory neurons",
#'                                                     "26"="Ciliated sensory neurons",
#'                                                     "27"="Oxygen sensory neurons",
#'                                                     "28"="Ciliated sensory neurons",
#'                                                     "29"="Unclassified neurons",
#'                                                     "30"="Socket cells",
#'                                                     "31"="Failed QC",
#'                                                     "32"="Pharyngeal gland",
#'                                                     "33"="Ciliated sensory neurons",
#'                                                     "34"="Ciliated sensory neurons",
#'                                                     "35"="Ciliated sensory neurons",
#'                                                     "36"="Failed QC",
#'                                                     "37"="Ciliated sensory neurons",
#'                                                     "38"="Pharyngeal muscle")
#'     neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
#'     pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn")
#'     pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#'     gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)
#'   }
#'
#' @export
find_gene_modules <- function(cds,
                          reduction_method = c("UMAP"),
                          max_components = 2,
                          umap.metric = "cosine",
                          umap.min_dist = 0.1,
                          umap.n_neighbors = 15L,
                          umap.fast_sgd = FALSE,
                          umap.nn_method = "annoy",
                          k = 20,
                          leiden_iter = 1,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = 0L,
                          cores=1,
                          verbose = FALSE,
                          preprocess_method = c('PCA', 'LSI'),
                          nn_control = list(),
                          ...) {
  method = 'leiden'

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=3,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=NULL,
                               k=k,
                               verbose=verbose)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  preprocess_method <- match.arg(preprocess_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(leiden_iter))
  ## TO DO what is resolution?
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))

  # preprocess_mat is gene_loading matrix. The gene_loadings were calculated only
  # for preprocess_method='PCA' in preprocess_cds() but I extend this to 'LSI' and
  # calculate gene_loadings here.
  preprocess_mat <- cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_v %*% diag(cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_sdev)

# Notes:
#   o  the beta vector is in cds@reduce_dim_aux[['Aligned']][['model']][['beta']]
#   o  cds@reduce_dim_aux[['Aligned']][['model']][['beta']] is npc x nfactor, which causes
#      preprocess_mat to have nfactor columns, often one column
#   o  I do not know how to adjust gene_loadings for batch effects
#      so this is disabled for now
#  if (!is.null(cds@reduce_dim_aux[['Aligned']][['model']][['beta']])){
#    preprocess_mat = preprocess_mat %*% (-cds@reduce_dim_aux[['Aligned']][['model']][['beta']])
#  }
  preprocess_mat <- preprocess_mat[intersect(rownames(cds), row.names(preprocess_mat)),]

  # uwot::umap uses a random number generator
  if( random_seed != 0L )
    set.seed( random_seed )

  umap_res = uwot::umap(as.matrix(preprocess_mat),
                        n_components = max_components,
                        metric = umap.metric,
                        min_dist = umap.min_dist,
                        n_neighbors = umap.n_neighbors,
                        fast_sgd = umap.fast_sgd,
                        n_threads=cores,
                        verbose=verbose,
                        nn_method= umap.nn_method,
                        ...)

  row.names(umap_res) <- row.names(preprocess_mat)
  if(ncol(umap_res) < 1) warning('bad loop: ncol(umap_res) < 1')
  colnames(umap_res) <- paste0('dim_', 1:ncol(umap_res))
  reduced_dim_res <- umap_res

  if(verbose)
    message("Running leiden clustering algorithm ...")

  cluster_result <- leiden_clustering(data=reduced_dim_res,
                                      pd=rowData(cds)[row.names(reduced_dim_res),,drop=FALSE],
                                      weight=weight,
                                      nn_index=NULL,
                                      k=k,
                                      nn_control=nn_control,
                                      num_iter=leiden_iter,
                                      resolution_parameter=resolution,
                                      random_seed=random_seed,
                                      verbose=verbose,
                                      ...)

  cluster_graph_res <- compute_partitions(cluster_result$g,
                                          cluster_result$optim_res,
                                          partition_qval, verbose)
  partitions <-
    igraph::components(cluster_graph_res$cluster_g)$membership[
      cluster_result$optim_res$membership]
  names(partitions) <- row.names(reduced_dim_res)
  partitions <- as.factor(partitions)

  gene_module_df <- tibble::tibble(id = row.names(preprocess_mat),
                                   module = factor(
                                     igraph::membership(cluster_result$optim_res)),
                                   supermodule = partitions)
  gene_module_df <- tibble::as_tibble(cbind(gene_module_df, umap_res))

  return(gene_module_df)
}





#' Aggregate columns within a sparse matrix.
#' @noRd
my.aggregate.Matrix = function (x, groupings = NULL, form = NULL, fun = "sum", ...)
{
  if (!methods::is(x, "Matrix"))
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  if (fun == "count")
    x <- x != 0
  groupings2 <- data.frame(A=as.factor(groupings))
  if (is.null(form))
    form <- stats::as.formula("~0+.")
  form <- stats::as.formula(form)
  mapping <- my.dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- Matrix::t(mapping) %*% x
  if (fun == "mean")
    result <- result/as.numeric(table(groupings)[rownames(result)])
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
                                                             groupings2$A))
  return(result)
}


# my.dMcast comes largely from the Matrix.utils::dMcast() function, which has
# been archived on CRAN because it appears to be unmaintained and it throws
# regression test errors.
#' Reshape data frame to a sparse matrix.
#' @noRd
my.dMcast <- function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  # values is a factor that multiplies the value called result, near the end of the function.
  # The multiplication method appears to be substantially more than a arithmetic multiplication.
  # The documentation for value.var says
  # @param value.var name of column that stores values to be aggregated numerics.
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]

  # all of the variables in the formula that relate to the data
  # terms returns a 'terms.object'
  # the terms.object help includes
  # Description: An object of class ‘terms’ holds information about
  #   a model. Usually the model was specified in terms of a ‘formula’
  #   and that formula was used to determine the terms object.
  # Value: The object itself is simply the formula supplied to the call
  # of ‘terms.formula’. The object has a number of attributes and they
  # are used to construct the model frame: factors, term.labels,
  # variables, intercept, order, response, offset, specials,
  # dataClasses, and predvars.
  alltms<-terms(formula,data=data)

  # get the name of the response variable that was givin in the
  # formula
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]

  # get the names of the explanatory variables
  tm<-attr(alltms,"term.labels")

  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))

  # get a character vector of all (variable) names in the
  # original formula
  allvars<-all.vars(alltms)

  # reshape the input data frame into a 'long' form that
  # has columns for the variables given in the original
  # formula
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)

  # make boolean vector describing whether or not the column is a factor (is this a correct description?)
  factors<-unlist(lapply(data,is.factor))

  # Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x)
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
        if(nlevels(x)!=length(na.omit(unique(x))))
          x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })

  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-Matrix::sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE)
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE)
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })

  # Make a dgCMatrix from the model matrix and response variables.
  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    # aggregate the matrix
    result<-my.aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}


#' Creates a matrix with aggregated expression values for arbitrary groups of
#' genes
#'
#' @param cds The cell_data_set on which this function operates
#' @param gene_group_df A dataframe in which the first column contains gene ids
#'   or short gene names and the second contains groups. If NULL, genes are not
#'   grouped.
#' @param cell_group_df A dataframe in which the first column contains cell ids
#'   and the second contains groups. If NULL, cells are not grouped.
#' @param norm_method How to transform gene expression values before
#'   aggregating them. If "log", a pseudocount is added. If "size_only", values
#'   are divided by cell size factors prior to aggregation.
#' @param pseudocount Value to add to expression prior to log transformation
#'   and aggregation.
#' @param scale_agg_values Whether to center and scale aggregated groups of
#'   genes.
#' @param max_agg_value If scale_agg_values is TRUE, the maximum value the
#'   resulting Z scores can take. Higher values are capped at this threshold.
#' @param min_agg_value If scale_agg_values is TRUE, the minimum value the
#'   resulting Z scores can take. Lower values are capped at this threshold.
#' @param exclude.na Logical indicating whether or not to exclude NA values
#'   from the aggregated matrix.
#' @param gene_agg_fun Function used for gene aggregation. This
#'   can be either sum or mean. Default is sum.
#' @param cell_agg_fun Function used for cell aggregation. Default is mean.
#'
#' @return A matrix of dimension NxM, where N is the number of gene groups and
#'   M is the number of cell groups.
#' @examples
#'   \dontrun{
#'      expression_matrix <- readRDS(system.file('extdata',
#'                                                'worm_l2/worm_l2_expression_matrix.rds',
#'                                                package='monocle3'))
#'      cell_metadata <- readRDS(system.file('extdata',
#'                               'worm_l2/worm_l2_coldata.rds',
#'                                package='monocle3'))
#'      gene_metadata <- readRDS(system.file('extdata',
#'                               'worm_l2/worm_l2_rowdata.rds',
#'                               package='monocle3'))
#'
#'      cds <- new_cell_data_set(expression_data=expression_matrix,
#'                               cell_metadata=cell_metadata,
#'                               gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds, num_dim = 100)
#'     cds <- reduce_dimension(cds)
#'     cds <- cluster_cells(cds, resolution=1e-5)
#'     colData(cds)$assigned_cell_type <- as.character(partitions(cds))
#'     colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
#'                                                     "1"="Germline",
#'                                                     "2"="Body wall muscle",
#'                                                     "3"="Unclassified neurons",
#'                                                     "4"="Vulval precursors",
#'                                                     "5"="Failed QC",
#'                                                     "6"="Seam cells",
#'                                                     "7"="Pharyngeal epithelia",
#'                                                     "8"="Coelomocytes",
#'                                                     "9"="Am/PH sheath cells",
#'                                                     "10"="Failed QC",
#'                                                     "11"="Touch receptor neurons",
#'                                                     "12"="Intestinal/rectal muscle",
#'                                                     "13"="Pharyngeal neurons",
#'                                                     "14"="NA",
#'                                                     "15"="flp-1(+) interneurons",
#'                                                     "16"="Canal associated neurons",
#'                                                     "17"="Ciliated sensory neurons",
#'                                                     "18"="Other interneurons",
#'                                                     "19"="Pharyngeal gland",
#'                                                     "20"="Failed QC",
#'                                                     "21"="Ciliated sensory neurons",
#'                                                     "22"="Oxygen sensory neurons",
#'                                                     "23"="Ciliated sensory neurons",
#'                                                     "24"="Ciliated sensory neurons",
#'                                                     "25"="Ciliated sensory neurons",
#'                                                     "26"="Ciliated sensory neurons",
#'                                                     "27"="Oxygen sensory neurons",
#'                                                     "28"="Ciliated sensory neurons",
#'                                                     "29"="Unclassified neurons",
#'                                                     "30"="Socket cells",
#'                                                     "31"="Failed QC",
#'                                                     "32"="Pharyngeal gland",
#'                                                     "33"="Ciliated sensory neurons",
#'                                                     "34"="Ciliated sensory neurons",
#'                                                     "35"="Ciliated sensory neurons",
#'                                                     "36"="Failed QC",
#'                                                     "37"="Ciliated sensory neurons",
#'                                                     "38"="Pharyngeal muscle")
#'     neurons_cds <- cds[,grepl("neurons", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
#'     pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn")
#'     pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#'     gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)
#'     cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)),
#'                                     cell_group=partitions(cds)[colnames(neurons_cds)])
#'     agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
#'   }
#'
#' @export
aggregate_gene_expression <- function(cds,
                                      gene_group_df = NULL,
                                      cell_group_df = NULL,
                                      norm_method=c("log", "binary",
                                                    "size_only"),
                                      pseudocount=1,
                                      scale_agg_values=TRUE,
                                      max_agg_value=3,
                                      min_agg_value=-3,
                                      exclude.na=TRUE,
                                      gene_agg_fun="sum",
                                      cell_agg_fun="mean"){
  if (is.null(gene_group_df) && is.null(cell_group_df))
    stop("one of either gene_group_df or cell_group_df must not be NULL.")
  agg_mat <- normalized_counts(cds, norm_method=norm_method,
                               pseudocount=pseudocount)
  if (is.null(gene_group_df) == FALSE){
    gene_group_df <- as.data.frame(gene_group_df)
    gene_group_df <- gene_group_df[gene_group_df[,1] %in%
                                     fData(cds)$gene_short_name |
                                     gene_group_df[,1] %in%
                                     row.names(fData(cds)),,drop=FALSE]

    # Convert gene short names to rownames if necessary. The more
    # straightforward single call to recode took much longer.
    # Thanks to Christopher Johnstone who posted this on github.
    short_name_mask <- gene_group_df[[1]] %in% fData(cds)$gene_short_name
    if (any(short_name_mask)) {
      geneids <- as.character(gene_group_df[[1]])
      geneids[short_name_mask] <- row.names(fData(cds))[match(
                  geneids[short_name_mask], fData(cds)$gene_short_name)]
      gene_group_df[[1]] <- geneids
    }

    # gene_group_df = gene_group_df[row.names(fData(cds)),]

    # FIXME: this should allow genes to be part of multiple groups. group_by
    # over the second column with a call to colSum should do it.
    gene_groups = unique(gene_group_df[,2])
    agg_gene_groups = lapply(gene_groups, function(gene_group){
      genes_in_group = unique(gene_group_df[gene_group_df[,2] == gene_group,1])
      gene_expr_mat = agg_mat[genes_in_group,]
      if (length(dn <- dim(gene_expr_mat)) < 2L)
        return(NA)
      if (gene_agg_fun == "mean"){
        res = Matrix::colMeans(agg_mat[genes_in_group,])
      }else if (gene_agg_fun == "sum"){
        res = Matrix::colSums(agg_mat[genes_in_group,])
      }
      return(res)
    })

    agg_mat_colnames = colnames(agg_mat)
    agg_mat = do.call(rbind, agg_gene_groups)
    row.names(agg_mat) = gene_groups
    agg_mat = agg_mat[is.na(agg_gene_groups) == FALSE,]
    colnames(agg_mat) = agg_mat_colnames
  }

  if (is.null(cell_group_df) == FALSE){

    cell_group_df <- as.data.frame(cell_group_df)
    cell_group_df <- cell_group_df[cell_group_df[,1] %in% row.names(pData(cds)),,
                                  drop=FALSE]
    agg_mat <- agg_mat[,cell_group_df[,1],drop=FALSE]
    agg_mat <- my.aggregate.Matrix(Matrix::t(agg_mat),
                                  as.factor(cell_group_df[,2]),
                                  fun=cell_agg_fun)
    agg_mat <- Matrix::t(agg_mat)
  }

  if (scale_agg_values){
    agg_mat_rownames = row.names(agg_mat)
    agg_mat_colnames = colnames(agg_mat)

    agg_mat = as.matrix(agg_mat)
    agg_mat <- t(scale(t(agg_mat)))
    agg_mat[is.nan(agg_mat)] = 0
    agg_mat[agg_mat < min_agg_value] <- min_agg_value
    agg_mat[agg_mat > max_agg_value] <- max_agg_value

    row.names(agg_mat) = agg_mat_rownames
    colnames(agg_mat) = agg_mat_colnames
  }

  if (exclude.na){
    agg_mat <- agg_mat[row.names(agg_mat) != "NA", colnames(agg_mat) != "NA",drop=FALSE]
  }
  return(agg_mat)
}
