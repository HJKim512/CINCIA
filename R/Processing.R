
#' Plotting gene/UMI counts for quality control
#'
#' Jitter plots and histograms of gene and UMI counts are created.
#' User can determine lower/upper thresholds for gene and UMI counts to filter out.
#'
#' @param rawObj A Seurat object of input gene expression matrix. This is not required if 'expMtx' is supplied.
#' @param expMtx Input gene expression count matrix. This is not required if 'rawObj' is supplied.
#' @param OutPath Path to the output directory under which subdirectories for saving the log files are created.
#' @param project.id A string which is used for 'orig.ident' in Seurat object and library name. Default to 'Project'.
#' @export

Lower_plot_Preprocess <- function( rawObj, expMtx, OutPath, project.id='Project'  ){
  require(ggplot2); require(Seurat)

  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  dir.create( paste0(OutPath, '1.Processing/raw'), recursive = TRUE )
  dir.create( paste0(OutPath, '1.Processing/stats'))
  dir.create( paste0(OutPath, '1.Processing/tmp'))

  if (missing(expMtx) & missing(rawObj)) stop('Either rawObj or expMtx should be supplied.')
  if (missing(expMtx)) {
    warning('\'expMtx\' was not supplied; Use rawObj instead..')

    cat(paste0(" Saving Expression matrix to the path : ", OutPath, '1.Processing/raw\n'))
    expMtx <- as.matrix(rawObj@assays$RNA@counts)
    expMtx.nonz <- expMtx[ rowSums(expMtx)!=0, ]
    write.table(expMtx.nonz, paste0(OutPath, "1.Processing/raw/rawData.ExpMtx.nonZero.txt"),  sep='\t', quote=FALSE)

    cat(" Saving label table to the path : ", OutPath, '\n\n')
    label <- data.frame(Cell_barcode = colnames(expMtx.nonz), library = rep(project.id, ncol(expMtx.nonz)), stringsAsFactors = FALSE)
    rownames(label) <- label$Cell_barcode
    write.table(label, paste0( OutPath, "1.Processing/raw/rawData.LabelTable.txt" ),  sep='\t', quote=FALSE, row.names = FALSE)

    cat(" Saving Seurat object to the path : ", OutPath, '\n\n')
    saveRDS(rawObj, paste0(OutPath, "1.Processing/tmp/rawObj.Rds"))
  }

  if (missing(project.id) & missing(rawObj)) warning('\'project.id\' was not supplied; using default project.id ' )
  if (missing(project.id) & missing(expMtx)) project.id <- levels(rawObj@meta.data$orig.ident)[1]

  if (missing(rawObj)) {
    warning('\'rawObj\' was not supplied; Use expMtx instead..')

    cat(" Saving Expression matrix to the path : ", OutPath, '\n')
    expMtx <- as.matrix(expMtx)
    expMtx.nonz <- expMtx[ rowSums(expMtx)!=0, ]
    #write.table(expMtx.nonz, paste0( OutPath, "1.Processing/raw/rawData.ExpMtx.nonZero.txt" ),  sep='\t', quote=FALSE)

    cat(" Saving label table to the path : ", OutPath, '\n\n')
    label <- data.frame(Cell_barcode = colnames(expMtx.nonz),
                        library = rep(project.id, ncol(expMtx.nonz)), stringsAsFactors = FALSE)
    rownames(label) <- label$Cell_barcode
    write.table(label, paste0( OutPath, "1.Processing/raw/rawData.LabelTable.txt" ),  sep='\t', quote=FALSE, row.names = FALSE)

    cat(" Creating and saving Seurat object to the path : ", OutPath, '\n\n')
    rawObj <- CreateSeuratObject(counts= expMtx.nonz, project=project.id)
    rawObj <- AddMetaData(object=rawObj, metadata=label)
    saveRDS(rawObj, paste0(OutPath, "1.Processing/tmp/rawObj.Rds"))
  }

  cat('Plotting jitter plots..\n')

  # UMI count
  plt1 <- ggplot(rawObj@meta.data, aes(orig.ident, nCount_RNA)) +
    geom_jitter(size = 0.25) + theme_bw() + xlab(label = '')+ ylab(label="UMI count" )+
    theme(text = element_text(size = 7), axis.text = element_text(size = 7),
          panel.grid = element_blank(), axis.text.x = element_text( vjust = .5))
  plt2 <- ggplot(rawObj@meta.data, aes(orig.ident, nCount_RNA)) +
    geom_jitter(size = 0.25) + ylim(0,2000)+ theme_bw() + xlab(label = '')+ ylab(label="UMI count" )+
    theme(text = element_text(size = 7), axis.text = element_text(size = 7),
          panel.grid = element_blank(), axis.text.x = element_text( vjust = .5))
  plt <- CombinePlots(list(plt1,plt2), ncol = 2)
  ggsave(paste0(OutPath, '1.Processing/stats/1.Jitter.UMIcnt.pdf'), units = 'cm', width = 10, height = 8)

  # Gene count
  plt1 <- ggplot(rawObj@meta.data, aes(orig.ident, nFeature_RNA)) +
    geom_jitter(size = 0.25) + theme_bw() + xlab(label = '') + ylab(label="Gene count" )+
    theme(text = element_text(size = 7), axis.text = element_text(size = 7),
          panel.grid = element_blank(), axis.text.x = element_text( vjust = .5))
  plt2 <- ggplot(rawObj@meta.data, aes(orig.ident, nFeature_RNA)) +
    geom_jitter(size = 0.25) + ylim(0, 2000) +
    theme_bw() + xlab(label = '')+ ylab(label="Gene count" )+
    theme(text = element_text(size = 7), axis.text = element_text(size = 7),
          panel.grid = element_blank(), axis.text.x = element_text( vjust = .5));
  plt <- CombinePlots(list(plt1,plt2), ncol = 2)
  ggsave(paste0(OutPath, '1.Processing/stats/1.Jitter.Genecnt.pdf'), units = 'cm', width = 10, height = 8)

  cat('Plotting Histograms..\n')

  ## Histograms
  # Gene count
  plt1 <- ggplot(rawObj@meta.data, aes(nFeature_RNA, col = orig.ident, fill= orig.ident)) +
    geom_histogram(binwidth = 1) + labs(x="Gene count", y="Cell count", title="Histograms of Gene count (all)") + theme_bw() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 10), plot.title = element_text(size = 10),
          panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + NoLegend()
  plt2 <- ggplot(rawObj@meta.data, aes(nFeature_RNA, col = orig.ident, fill= orig.ident)) +
    geom_histogram(binwidth = 1) + labs(x="Gene count", y="Cell count", title="Histogram of Gene count (low)") + theme_bw() +
    xlim(0, 2000) + geom_vline(xintercept = c(100, 200, 300, 400, 500), col = 'black') +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 10), plot.title = element_text(size = 10),
          panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + NoLegend()
  plt <- CombinePlots(list(plt1,plt2), ncol = 1);plt
  ggsave(paste0(OutPath, '1.Processing/stats/1.Histogram.Genecnt.pdf'), units = 'cm', width = 30, height = 20)

  # UMI count
  plt1 <- ggplot(rawObj@meta.data, aes(nCount_RNA, col = "#2a9d8f", fill= orig.ident)) +
    geom_histogram(binwidth = 1, color="#2a9d8f") + labs(x="UMI count", y="Cell count", title="Histograms of UMI count (all)") +
    theme_bw() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 10),panel.grid = element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), plot.title = element_text(size = 10))+ NoLegend()

  plt2 <- ggplot(rawObj@meta.data, aes(nCount_RNA, col = "#2a9d8f", fill= orig.ident)) +
    geom_histogram(binwidth = 1, color="#2a9d8f") + labs(x="UMI count", y="Cell count",  title="Histogram of UMI count (low)") + theme_bw() +
    xlim(0, 2000) + geom_vline(xintercept = c(500, 600, 700, 800, 900, 1000), col = 'black') +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 10), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), plot.title = element_text(size = 10) ) + NoLegend()
  plt <- CombinePlots(list(plt1,plt2), ncol = 1);plt
  ggsave(paste0(OutPath, '1.Processing/stats/1.Histogram.UMIcnt.pdf'), units = 'cm', width = 30, height = 20)

  cat('Plots were saved. Determine thresholds for filtering out gene/UMI counts. \n\n')
}




#' Filtering out cells based on predetermined lower/upper thresholds for gene/UMI counts.
#'
#' Cells with too low gene/UMI counts are filtered out based on user-defined thresholds (thres.gene, thres.umi) in previous step.
#' For upper gene count filtering, either thres.UpGene or thres.prob should be provided.
#' If user put in thres.prob (in range 0 to 1), parameters of Gaussian mixture model (number of components = 3)
#' fitted gene counts are estimated, and cells with higher posterior probability of belonging to component-3 than thres.prob (default=0.7)
#' are filtered out, being regarded as putative multiplets.
#'
#' @param rawObj A Seurat object created with input expression matrix.
#' @param thres.gene A numberic value indicating lower gene count threshold.
#' @param thres.umi A numberic value indicating lower UMI count threshold.
#' @param thres.UpGene A numberic value indicating upper gene count threshold. Either thres.UpGene or thres.prob should be supplied.
#' @param thres.prob Threshold for posterior probability of cells belonging to component-3 in the gene count mixture model.
#' Cells with their posterior probability higher than this value (Defaults to 0.7) are regarded as putative multiplets and filtered out.
#' Either thres.UpGene or thres.prob should be supplied.
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' @export

Upper_Preprocess <- function( rawObj, thres.gene, thres.umi, thres.UpGene, thres.prob=0.9, OutPath  ){

  if (missing(thres.gene)) stop('\'thres.gene\' was not supplied.')
  if (missing(thres.umi)) stop('\'thres.umi\' was not supplied.')
  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  require(mixtools); require(ggplot2); require(Seurat); require(dplyr)

  if (!missing(thres.UpGene)){
    Obj.flt <- subset(rawObj, subset = nFeature_RNA > thres.gene & nCount_RNA > thres.umi & nFeature_RNA < thres.UpGene)
    cat("\nSaving lower/upper gene/UMI-filtered object ..\n")
    saveRDS(Obj.flt, paste0(OutPath, "1.Processing/tmp/rawObj.AllFlt.Rds"))
    return(Obj.flt)
  }
  else{
    warning('\'thres.UpGene\' was not supplied; Using \'thres.prob\' instead ..')

    Obj.flt <- subset(rawObj, subset = nFeature_RNA >= thres.gene & nCount_RNA >= thres.umi )
    cat("\nSaving lower gene-filtered object ..\n")
    saveRDS(Obj.flt, paste0(OutPath, "1.Processing/tmp/rawObj.lowerFlt.Rds"))

    # Upper gene filtering
      # define mix_comps()
    mix_comps <- function(x, mu, sigma, lam){ lam * dnorm(x, mu, sigma) }
    mixmdl <- normalmixEM(Obj.flt$nFeature_RNA, k = 3)

    plt <- data.frame(x = mixmdl$x) %>%
      ggplot() + geom_histogram(aes(x, ..density..), binwidth = 1, colour = "#B6CEC7", fill = "white") +
      stat_function(geom = "line", fun = mix_comps, args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "#ED5564", lwd = 1.2) +
      stat_function(geom = "line", fun = mix_comps, args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "#FFCE54", lwd = 1.2) +
      stat_function(geom = "line", fun = mix_comps, args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
                    colour = "#019875", lwd = 1.2) +
      theme_bw()+ labs(title = "Gene count of the mixture model", y="Density", x="Gene count")+
      theme(plot.title = element_text(hjust=0.5, size=22), axis.title = element_text(size=16), axis.text = element_text(size=15))
    ggsave( paste0(OutPath, "1.Processing/stats/1.Histogram.GMM.pdf"), units = 'cm', height=18, width=42)

    post.df <- data.frame(row.names = rownames(Obj.flt@meta.data), as.data.frame(cbind(x_1 = mixmdl$x, mixmdl$posterior)))

    cat( " Example : p(x|comp.3) > ", thres.prob , '\n\n')
    ex <- post.df[post.df$comp.3 >= thres.prob, ]
    print( ex[order(ex$x),][1:20,])

    cut <- rownames(post.df[post.df$comp.3 > thres.prob,])

    Obj.upflt <- subset(Obj.flt, cells= setdiff(colnames(Obj.flt), cut))
    cat("\nSaving upper gene-filtered object ..\n")
    saveRDS(Obj.upflt, paste0(OutPath, "1.Processing/tmp/rawObj.AllFlt.Rds") )

    return(Obj.upflt)
  }
}


#' Calculate mitochondrial gene expression and plotting the distribution.
#'
#' Percentage of mitochondrial gene expression in each cell is calculated and plotted.
#'
#' @param Obj.flt A Seurat object where cells with too low gene/UMI counts or high gene counts are filtered out.
#' @param prefix Prefix of mitochondrial genes of the species of your data. (not required but recommended)
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' @export

mtRNA_plot_Preprocess <- function( Obj.flt, prefix, OutPath ){

  if (missing(Obj.flt)) stop('\'Obj.flt\' was not supplied.')

  if (missing(prefix)){
    warning('\'prefix\' was not supplied. Searching for the prefix of genes encoding mitochondrial RNAs..')
    if ( length( grep(rownames(Obj.flt), value=TRUE, pattern='^MT-')>0 )) prefix <- 'MT-'
    else if ( length(grep(rownames(Obj.flt), value=TRUE, pattern='^mt-'))>0 ) prefix <- 'mt-'
    else prefix <- 'mt:'
    cat("The prefix of gene is ", prefix, '\n\n')
  }

  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  require(ggplot2);  require(Seurat)

  Obj.flt[["percent.mt"]] <- PercentageFeatureSet(object = Obj.flt, pattern = prefix)

  plt1 <- ggplot(Obj.flt@meta.data, aes(percent.mt, col = orig.ident)) + geom_density() + NoLegend()+
    theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7),
                       panel.grid = element_blank(), legend.position = 'none')
  plt2 <- ggplot(Obj.flt@meta.data, aes(percent.mt, fill = orig.ident, color='#6D8700')) +
    geom_histogram(binwidth = 1) + NoLegend()+
    theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7),
                       panel.grid = element_blank(), legend.position = 'none')
  plt <- CombinePlots(list(plt1,plt2), ncol = 1)
  ggsave(paste0(OutPath, '1.Processing/stats/1.Plots.percent.mt.pdf'), units = 'cm', width = 16, height = 10)
  cat("\nSaving mitoGene plot ..\n")
  return(Obj.flt)
}


#' Normalize the data, find highly variable genes, scale the data and perform PCA
#'
#' Normalization, Finding variable genes, scaling and PCA are performed.
#'
#' @param Obj.flt A Seurat object where cells with too low gene/UMI counts or high gene counts are filtered out.
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' @param num.dim Number of principal components to calculate. Defaults to the largest available limit, 50.
#' @export
NormToPCA_Process <- function( Obj.flt, OutPath, num.dim=50 ){
  if (missing(Obj.flt)) stop('\'Obj.flt\' was not supplied.')

  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  if (missing(num.dim)) warning('\'num.dim\' was not supplied; Using default 50 ..')

  require(ggplot2); require(Seurat)

  Obj.norm <- NormalizeData(object = Obj.flt, normalization.method = "LogNormalize", scale.factor = 10000)
  Obj.norm <- FindVariableFeatures(object = Obj.norm, selection.method = "vst")

  top10 <- head(x = VariableFeatures(object = Obj.norm), 10)
  plot1 <- VariableFeaturePlot(object = Obj.norm)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(paste0(OutPath, '1.Processing/stats/2.FindVariableFeatures.mitoCut.pdf'), units = 'cm', width = 17, height = 10)

  Obj.norm <- ScaleData(object = Obj.norm, vars.to.regress = c("percent.mt", "nCount_RNA"))
  Obj.norm <- RunPCA(object = Obj.norm, features = VariableFeatures(Obj.norm))
  saveRDS(Obj.norm, paste0(OutPath, '1.Processing/tmp/rawObj.norm.PCA.Rds'))

  Obj.norm <- JackStraw(object = Obj.norm, num.replicate = 50, dims = num.dim)
  Obj.norm <- ScoreJackStraw(object = Obj.norm, dims = 1:num.dim)
  JackStrawPlot(Obj.norm, dims = 1:num.dim)
  ggsave(paste0(OutPath,'1.Processing/stats/2.pca.JackStrawPlot.mitoCut.pdf'), units = 'cm', width = 40, height = 12)
  saveRDS(Obj.norm, paste0(OutPath, '1.Processing/tmp/rawObj.norm.PCA.Rds'))

  return(Obj.norm)
}




#' Running UMAP, Louvain algorithm for clustering and Finding differentially expressed genes (DEGs)
#'
#' UMAP, clustering and DE analysis are done to the input data sequentially, a Seurat object on which PCA was previously performed .
#'
#' @param Obj.norm A Seurat object where analysis processes are all performed in the previous step 'NormtoPCA_Process()'.
#' @param nPC Number of principal components to use. Default to 20.
#' @param genes A vector containing genes to use when plotting features.
#' @param seed Random seed number to use for generating UMAP plots. Default to '777778'
#' @param res Value of the resolution parameter which is proportional to the number of communities. Default to 0.8
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' @export
#'
UMAPtoDEG_Process <- function( Obj.norm, nPC=20, genes, seed=777778, res=0.8, OutPath ){
  if (missing(Obj.norm)) stop('\'Obj.norm\' was not supplied.')
  if (missing(nPC)) warning('\'nPC\' was not supplied; Using default 20')
  if (missing(genes)) stop('\'genes\' was not supplied.')

  if (missing(seed)){
    warning('\'seed\' was not supplied; Using default 777777')
  }

  if (missing(res)){
    warning('\'res\' was not supplied; Using default 0.8')
  }

  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  dir.create( paste0(OutPath, '1.Processing/deg') )

  require(MAST); require(ggplot2); require(Seurat); require(dplyr)

  # UMAP
  cat("Running UMAP ...\n\n")
  Obj.norm <- RunUMAP(Obj.norm, dims=1:nPC, reduction = "pca", reduction.key='UMAP',
                      n.components=2, min.dist=0.2, seed.use = seed)
  plt1 <- FeaturePlot(Obj.norm, features = genes, cols = c("grey","red"), reduction = "umap")
  ggsave(paste0(OutPath, "1.Processing/stats/2.UMAP.markers.seed.", seed, ".pdf"), units = 'cm', width = 48, height = 48)

  # Clustering
  cat("Louvain Clustering ...\n\n")
  Obj.norm <- ScaleData(object = Obj.norm)
  Obj.norm <- FindNeighbors(object = Obj.norm, dims = 1:nPC)
  Obj.norm <- FindClusters(object = Obj.norm, resolution = res)

  plt2 <- DimPlot(Obj.norm, reduction = "umap", pt.size = .05, label = T, label.size = 2)+
    theme(axis.title = element_text(size=7), axis.text = element_text(size=7),
          legend.text = element_text(size=9), legend.key.width = unit(0.2,'cm'), legend.key.height=unit(0.2,'cm') )
  ggsave(paste0(OutPath, "1.Processing/stats/2.Cluster.min_dist.0.2.res.", res,".pdf"), units = 'cm', width = 13, height = 10)
  saveRDS(Obj.norm, paste0(OutPath, "1.Processing/tmp/rawObj.Cluster.Rds"))

  # Find DEGs
  cat("Finding DEGs ...\n\n")
  Obj.norm.markers <- FindAllMarkers(object = Obj.norm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
  top10 <- Obj.norm.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC); top10
  dehm <- DoHeatmap(object = Obj.norm, features = top10$gene, angle = 90, size = 1, raster = T, draw.lines = T, label=T)
  ggsave(paste0(OutPath, "1.Processing/deg/1.Heatmap.markers.res.", res,".pdf"), units = 'cm', width = 50, height = 40)
  write.table(Obj.norm.markers, paste0(OutPath, "1.Processing/deg/2.Table.markers.res", res,".txt"), sep = '\t', quote = F, col.names = T, row.names = F)

  return(Obj.norm)
}



#' Running DoubletFinder (C. S. McGinnis et.al 2019)
#'
#' @param object A fully-processed (from QC to cell type annotation) Seurat object.
#' @param npc Number of principal components to use.
#' @param nLib The number of libraries (samples) of data to run. Default to 1
#' @param nRun The number of running iteration with different random seed. Default to 1
#' @param rseed Initial random seed to use for running DoubletFinder. Each of nRun iterations uses rseed multiplied by iteration number.
#' @param est.rates Vector containing estimated rate of doublets (multiplets) of the data, which is in range (0,1).
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' @export
Run_DoubletFinder <- function( object, npc , nLib=1, nRun=1, rseed=2022, est.rates, OutPath){

  if (missing(object)) stop('\'object\' was not supplied')
  if (missing(npc)) warning('\'nPC\' was not supplied; Using default 20')
  if (missing(nLib)) warning('\'nLib\' was not supplied; Using default 1')
  if (missing(nRun)) warning('\'nRun\' was not supplied; Using default 1')
  if (missing(rseed)) warning('\'rseed\' was not supplied; Using default 2022')
  if (missing(est.rates)) stop('Vector \'est.rates\' was not supplied')

  if (missing(OutPath)) stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  require(DoubletFinder); require(Seurat); require(dplyr)

  resFiles <- list()

  cat("1  Pre-process the input object...\n\n")
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures( object, selection.method = "vst", nfeatures = 2000)
  object <- ScaleData(object = object )
  object <- RunPCA(object, features = VariableFeatures(object), npcs = npc)
  object <- RunUMAP(object, dims = 1:npc, min.dist=0.2 )

  cat("2  Running DoubletFinder ...\n\n")
  for (i in 1:nLib) {

    # Basic inputs
    cat(paste('Current iterations : sample ',i,'\n',sep=''))

    homotypic.prop <- modelHomotypic(object@meta.data$celltype)
    nExp_poi <- round( est.rates[i] * nrow(object@meta.data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    # multiple times
    for (j in 1:nRun){
      cat(paste('iteration ', j, '\n', sep=''))
      seed <- rep(rseed*j)
      set.seed(seed)

      sweep.res.list <- paramSweep_v3(object, PCs = 1:npc, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- as.data.frame(find.pK(sweep.stats))
      select_pK <- as.numeric(as.character(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),'pK']))

      object <- doubletFinder_v3(object, PCs = 1:npc, pN = 0.25, pK = select_pK, nExp = nExp_poi,
                                 reuse.pANN = FALSE, sct = FALSE)
      reuse <- grep(colnames(object@meta.data), value = TRUE, pattern = paste("pANN_0.25", select_pK, nExp_poi,sep="_"))

      object <- doubletFinder_v3(object, PCs = 1:npc, pN = 0.25, pK = select_pK, nExp = nExp_poi.adj, reuse.pANN = reuse, sct = FALSE)
      class <- grep(colnames(object@meta.data), value=TRUE, pattern = paste('DF.classifications_0.25',select_pK,nExp_poi.adj,sep='_'))

      if(j==1){  res = data.frame( row.names = rownames(object@meta.data), object@meta.data[, class]) }
      else{ res = cbind(res, object@meta.data[, class]) }
    }
    colnames(res) <- rep(paste("iter", 1:nRun, sep="_"))
    lib_name <- paste("Library",i,sep="_")
    resFiles[[lib_name]] <- res
  }

  # generate singlet object
  if (nLib==1 & nRun==1){
    object$DoubletFinder <- resFiles$Library_1$iter_1
    object.sgl <- subset(object, DoubletFinder=="Singlet")

    cat("\nSaving object with class labels predicted from DoubletFinder..\n")
    saveRDS(object, paste0(OutPath, "1.Processing/tmp/rawObj.DFresult.Rds") )
    saveRDS(object.sgl, paste0(OutPath, "1.Processing/tmp/Reference.Set.Rds") )
    return(object.sgl)
  }
  else return(resFiles)
}


