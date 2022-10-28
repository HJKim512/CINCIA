
#' Assigning two representative cell types to each doublet
#'
#' For each doublet, two cell types of top 2 high prediction scores are assigned.
#'
#' @param predicted.df A data frame of cell type prediction scores for each putative doublet (cell barcode X cell types),
#' which is resulted from functions 'FindTransferAnchors' and 'TransferData'.
#' @export
FindTwoCelltypes <- function( predicted.df ){

  ct.pool <- rownames(predicted.df) # All the cell types in the data
  barcodes <- colnames(predicted.df) # Cell barcodes of doublets

  result <- data.frame( row.names = c(), top.ct1 = c(), top.ct2 = c() )

  for (i in 1:ncol(predicted.df)) {
      proportions <- predicted.df[[i]]
      tops <- proportions[ order(proportions, decreasing = TRUE) ][1:2]  # cell type proportions of top2 cell types
      ct <- ct.pool[ c( which(proportions==tops[1]), which(proportions==tops[2]) )]
      df <- data.frame(row.names=barcodes[i], top.ct1 = ct[1], top.ct2 = ct[2])
      result <- rbind(result, df)
  }
  return(result)
}



#' Cell type Assignment (Doublet deconvolution)
#'
#' Transferring cell types from reference set to doublet set using functions 'FindTransferAnchors' and 'TransferData' in Seurat pipeline.
#'
#' @param nPC The number of PCs to use for label transferring. Default to 20
#' @param refObj Reference set; A fully processed (from QC to annotation) Seurat object
#' @param rawObj A Seurat object of raw input expression matrix
#' @param NBres Returned object from function 'DefineDoublets', where putative doublets called from doublet meta-caller and
#' their information are included
#' @param thres Threshold value in range (0, 1), which is used to filter maximal cell type prediction score. Default to 0.98
#' @param OutPath Path to the output directory under which log files of defining doublet set are saved.
#' (Path used as the parameter PATH_out in previous step of Define_doublets().)
#' @param col.label Name of the column containing cell type label in the reference set
#' @export
#'
Assign_Celltypes <- function( refObj, rawObj, NBres, col.label, nPC=20, thres=0.98, OutPath ){

  if (missing(nPC)) warning('\'nPC\' was not supplied, using default 20')
  if (missing(refObj)) stop('\'refObj\' was not supplied')
  if (missing(rawObj)) stop('\'rawObj\' was not supplied')
  if (missing(NBres)) stop('\'NBres\' was not supplied')
  if (missing(thres)) warning('\'thres\' was not supplied, using default 0.98')
  if (missing(col.label)) stop('\'col.label\' was not supplied')
  if (missing(OutPath)) stop('\'OutPath\' was not supplied')

  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  rows <- c(); cell_1 <- c(); cell_2 <- c(); result <- c()
  require(Seurat); require(reshape2)

  # 1. Label transfer to the doublet object
  cat(" 1 Transfer cell type information from reference set to doublets ...\n\n")

  # Generate doublet set and process
  doubletset <- subset(rawObj, cells=rownames(NBres[NBres$predict=='doublet',]))
  doubletset <- NormalizeData(doubletset)
  doubletset <- FindVariableFeatures(doubletset)
  doubletset <- ScaleData(doubletset)

  anchors <- FindTransferAnchors(reference = refObj, query = doubletset,
                                 reference.reduction = 'pca', scale = FALSE, dims=1:nPC)
  prediction <- TransferData(anchorset = anchors, refdata = refObj[[col.label]][,col.label],
                             reference = refObj, n.trees = 100, dims = 1:nPC)
  # For each cell type, extract prediction score only
  prediction.df <- as.data.frame(t(prediction[,2:ncol(prediction)]))

  # Modify cell type names (remove prefix 'predicted.score.')
  for (i in 1:length(rownames(prediction.df))) {
    name <- strsplit( rownames(prediction.df)[i], split='score.' )[[1]][2]
    rows[i] <- name
  }
  rownames(prediction.df) <- rows

  prediction.df <- as.data.frame(t(prediction.df))
  cat("     The number of detected doublets : ", nrow(prediction.df))

  # Modify the cell type names (to be identical to original name)
  colnames(prediction.df) <- gsub('\\.', " ", colnames(prediction.df))
  prediction.df.nomax <- prediction.df[, 1:ncol(prediction.df)-1]  # exclude 'max' score column

  # Assign top 2 cell types to the doublets
  cat("\n\n 2 Assign cell types of the two highest scores to doublets ...\n")
  input <- as.data.frame(t(prediction.df.nomax))
  composition <- FindTwoCelltypes( input[1:(nrow(input)-1),] )

  # filtering maximum prediction score
  cat("\n\n 3 Filter doublets whose maximum prediction score is higher than ", thres, " ...\n")
  prediction.df.flt <- prediction.df[as.numeric(prediction.df$max) < thres, ]
  composition.flt <- composition[rownames(prediction.df.flt),]
  cat("\n     The number of doublets left : ", nrow(prediction.df.flt))

  # Count combination
  cat("\n\n 4 Count doublets belonging to each cell type combination ...\n")
  for (i in 1:nrow(composition.flt)) {
      ct.pool <- c( composition.flt[i,'top.ct1'], composition.flt[i,'top.ct2'] )
      ct.pool <- sort(ct.pool)
      composition.flt[i,'comb'] <- paste(ct.pool[1],ct.pool[2],sep='-')
  }

  # combination table
  DblClass <- as.data.frame(table(composition.flt$comb))
  rownames(DblClass) <- DblClass$Var1

  for (i in 1:nrow(DblClass)) {
      tmp <- strsplit(as.character(DblClass[i,'Var1']), '-')[[1]]
      cell_1[i] <- tmp[1]; cell_2[i] <- tmp[2]
  }
  DblClass$cell_1 <- cell_1;  DblClass$cell_2 <- cell_2
  DblClass <- DblClass[,2:4] # remove 'Var'

  cat("\n\n Finished \n")
  result <- list( 'all.predict'=prediction.df,
                  'all.composition'=composition,
                  'flt.predict'=prediction.df.flt,
                  'flt.composition'=composition.flt,
                  'combination'=DblClass )

  saveRDS(result, paste0(OutPath, "2.Doublet/3.CelltypeAssignment.Allresult.Rds"))
  return(result)
}



#' Hypergeometric test
#'
#' @param data A data frame containing information of doublets including potential interactions
#' @export
Hypertest <- function(data){
  result <- c()

  for (i in 1:nrow(data)) {
    overlap <- as.numeric(data[i,'Freq'])
    success_bg <- sum(data[data$cell_1==data[i, 'cell_2'] ,'Freq'])
    fail_bg <- sum(data$Freq)-success_bg
    sample <- sum(data[data$cell_1==data[i, 'cell_1'] | data$cell_2==data[i, 'cell_1'],'Freq'])

    phyper(overlap-1, success_bg, fail_bg, sample, lower.tail = FALSE)
    result[rownames(data[i,])] <- phyper(overlap-1, success_bg, fail_bg, sample, lower.tail = FALSE)
  }
  return(result)
}




#' Finding significant interactions
#'
#' Fraction of each cell type in the data are calculated first to estimate the expected number of doublets.
#' Then, observed number of doublets are compared to the expected number for each cell type combination (interaction).
#' P-values are calculated by hypergeometric test (equivalent to one-sided Fisher's exact test) and adjusted through Bonferroni correction.
#'
#' @param DeconvObj Result object from function 'Assign_Celltypes()'; a list containing doublet deconvolution result
#' @param refObj Reference set; A fully processed (from QC to annotation) Seurat object
#' @param col.label Name of the column containing cell type label in the reference set
#' @param OutPath Path to the output directory under which log files of defining doublet set are saved.
#' (Path used as the parameter PATH_out in previous step of Define_doublets().)
#' @param lower.frac The minimal fraction (prediction score calculated from Cell type Assignment) of a cell-type
#' to be considered as a constituent cell of a droplet (multiplet).
#' @export
#' @return A list containing cell type fraction of the data and information of all detected doublets including potential interactions
FindInteraction <- function( DeconvObj, refObj, col.label, lower.frac=0.1, OutPath ){

  if (missing(refObj)) stop('\'refObj\' was not supplied')
  if (missing(DeconvObj)) stop('\'DeconvObj\' was not supplied')
  if (missing(col.label)) stop('\'col.label\' was not supplied')
  if (missing(lower.frac)) warning('\'lower.frac\' was not supplied, using default 0.1 ')
  if (missing(OutPath)) stop('\'OutPath\' was not supplied')

  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')


  result <- c(); nMtl <- c(); nExp <- c()

  # Calculate cell type fraction in the data
  # calculate the number of cells associated in partial multiplets (using label transfer result)
  pred.scores <- DeconvObj[['flt.predict']]
  pred.scores <- pred.scores[, setdiff(colnames(pred.scores), 'max')]

  count <- as.data.frame(t(pred.scores)) %>%
    apply( ., MARGIN = 1, FUN = function(x){ ifelse(x>=lower.frac, 1, 0) } )
  sum.count <- colSums(count)

  # Combine with Cell type fractions in singlets (reference set)
  CelltypeProps <- data.frame( table(refObj[[col.label]]), row.names = 1)
  CelltypeProps <- data.frame( row.names = rownames(CelltypeProps),
                               'singlet'=CelltypeProps$Freq,
                               'doublet' = sum.count[rownames(CelltypeProps)])
  CelltypeProps$total <- CelltypeProps$singlet + CelltypeProps$doublet
  CelltypeProps <- data.frame(CelltypeProps, 'celltype' = rownames(CelltypeProps),
                              'fraction' = CelltypeProps$total/sum(CelltypeProps$total))

  # Calculate expected number of doublets and compare
  DblClass <- DeconvObj[['combination']]

  # find the number of doublets including cell_1
  for (i in 1:nrow(DblClass)){
      ct <- DblClass[i,'cell_1']
      nMtl[i] <- sum( DblClass[DblClass$cell_1==ct | DblClass$cell_2==ct,'Freq'])
      nExp[i] <- nMtl[i] * CelltypeProps[DblClass[i,'cell_2'], 'fraction']
  }
  DblClass$nMultiplet <- nMtl
  pseudocnt <- 1
  DblClass$ExpectedNum <- round(nExp) + pseudocnt
  DblClass$ratio <- DblClass$Freq/DblClass$ExpectedNum
  DblClass$pval <- Hypertest(DblClass)
  DblClass$BonF <- p.adjust(DblClass$pval, 'bonferroni')

  options(digits=8)

  result <- list( 'celltype.fraction'=CelltypeProps,
                  'combination.final'=DblClass )
  saveRDS(result, paste0(OutPath, "2.Doublet/3.Inference.Allresult.Rds"))

  signif <- result$combination.final[ (result$combination.final$BonF<0.05 &
                                       result$combination.final$ratio>1.0 ),]
  saveRDS(signif, paste0(OutPath, "2.Doublet/3.Inference.sig.Rds"))

  return(signif)
}











