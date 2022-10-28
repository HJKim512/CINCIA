
#' Heterotypic doublet simulation using reference set (singlets)
#'
#' Artificial heterotypic doublets are simulated using cells of reference Set, and combined with singlets together
#' into a gene expression matrix ( genes X cell barcodes ) and label table containing both singlets and doublets is generated.
#'
#' @param refObj Reference Set. A fully-processed Seurat object, containing putative singlets only.
#' @param nSimul The number of simulation runs. The number of simulated expression matrices as many as nSimul are created. Default to 100
#' @param nDbl The number of artificial doublets to generate for each simulation, which is recommended to be determined based on the number of cells
#' and multiplet rates of scRNA-seq method (ex. 10X Genomics Chromium) used. For each simulated dataset, artificial doublets within plus or minus 50 doublets
#' from nDbl are generated.
#' @param col.ref A Character string, name of the column which contains cell type label
#' @param OutPath Path to the output directory, where its subdirectories contain log files generated in previous steps.
#' New folders named 'simul_#'(total nSimul folders) are created and each simulated dataset is saved in there.
#' @export
#'
Simulate_Hetero_Doublets <- function( refObj, nSimul=100, nDbl, col.ref, OutPath ){
  all.result <- list()
  require(data.table)
  require(dplyr)

  if (missing(refObj))
      stop('\'refObj\' was not supplied : A fully-processed Seurat object containing singlets is required.')
  if (missing(nDbl))
      stop('\'nDbl\' was not supplied')

  if (missing(OutPath))
      stop('\'OutPath\' was not supplied')
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')

  if (missing(nSimul)) {
      warning('\'nSimul\' was not supplied, using default 100')
      nSimul = 100
  }
  if (missing(col.ref)) stop('\'col.ref\' was not supplied')

  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  dir.create( paste0(OutPath, '2.Simul/Input'), recursive=TRUE )
  dir.create( paste0(OutPath, '2.Simul/Output'), recursive=TRUE )

  expMtx <- refObj@assays$RNA@counts
  celltype.pool <- data.frame( row.names=colnames(refObj), refObj[[col.ref]] )
  #print(head(celltype.pool))
  nDbls <- sample( seq(nDbl-50, nDbl+50, by=2), nSimul, replace = TRUE ) # vector

  # 1. generate doublet sets and its reference

  for (i in 1:nSimul) {
    cat('# Simulation set ', i, '\n')
    cells.pool <- colnames(refObj)  # repeatedly update per every simulation set
    setname <- paste('simul', i, sep='_')
    bc.1=c(); bc.2=c(); ct.1=c(); ct.2=c()

    # randomly generate heterotypic doublets
    for (j in 1:nDbls[i]) {

      # sample two cells from two different cell types
      bc.1[j] <- base::sample(cells.pool, 1)
      ct.1[j] <- as.character(celltype.pool[bc.1[j], col.ref])

      subpool <- intersect( rownames(subset(celltype.pool, col.ref != ct.1[j])),  cells.pool )
      bc.2[j] <- base::sample(subpool, 1)
      ct.2[j] <- as.character(celltype.pool[bc.2[j], col.ref])

      # remove sampled cells from cells.pool
      cells.pool <- setdiff(cells.pool, c(bc.1[j], bc.2[j]))
    }

    res <- data.frame( row.names = paste('dbl', 1:nDbls[i], sep='_'), bc.1, bc.2, ct.1, ct.2 )
    all.result[[setname]] <- res

    # 2. Merge putative singlets and simulated doublets

    art.dbls <- as.data.frame( expMtx[, res$bc.1] + expMtx[, res$bc.2])
    colnames(art.dbls) <- rownames(res)
    simulMtx <- cbind(expMtx, art.dbls)
    simulMtx.nonZ <- simulMtx[ rowSums(simulMtx)!=0, ]

      # Label Table
    label <- data.frame(row.names = colnames(simulMtx.nonZ),
                        Cell_barcode = colnames(simulMtx.nonZ),
                        library = rep(as.character(refObj$orig.ident[1]), ncol(simulMtx.nonZ)), stringsAsFactors = FALSE)

      # save simulated expression matrix and label table as text files
    path <- paste0(OutPath, '2.Simul/Input/simul_', i )
    dir.create(path)

    cat('   - Saving expression matrix . . .', '\n')
    # 22.08.26 modified - fwrite
    fwrite(simulMtx.nonZ, paste( path, "/simulDbls.ExpMtx.txt", sep=''), sep='\t', quote = FALSE, row.names = T)
    #write.table(simulMtx.nonZ, paste( path, "/simulDbls.ExpMtx.txt", sep=''),  sep='\t', quote=FALSE)

    cat('   - Saving Label table . . .', '\n')
    # 22.08.26 modified - fst package
    fwrite(label, paste(path, "/simulDbls.LabelTable.txt", sep=''), sep='\t', quote = FALSE, row.names=T)
    #write.table(label, paste( path, "/simulDbls.LabelTable.txt", sep=''),  sep='\t', quote=FALSE, row.names = FALSE)
  }
  saveRDS(all.result, paste0(OutPath, '2.Simul/0.SimulDoubletInfo.Rds'))
}


#' Doublet scoring through Scrublet, a python-based algorithm (Wolock, et al. 2019)
#'
#' Doublet score (probability) of each cell in the input data or simulated data is yielded by Scrublet.
#'
#' @param PyPath Path to the Python executable in your computing environment (Python environment path).
#' Make sure the package 'Scrublet' is already installed.
#' @param InPath Path to the input directory in which input count matrix (for type='raw') or folders named 'simul_#' which contain
#' count matrices of doublet-simulated datasets (for type='simul') exist.
#' @param OutPath Path to the output directory in which the results files are saved.
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param nSet The number of dataset to run.
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param est.rate Estimated doublet rates in your dataset.
#' @export
DoubletScoring_Scrublet <- function( PyPath, InPath, OutPath, type, nSet, nRun=10, est.rate=NULL){

  if (missing(PyPath)) stop("\'PyPath\', the path to the python you want to use was not supplied")
  if (missing(InPath)) stop("\'InPath\' was not supplied")
  if (grepl("/$", InPath) == FALSE) InPath <- paste0(InPath, '/')
  if (missing(OutPath)) stop("\'OutPath\' was not supplied")
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  if (missing(type)) stop("\'type\' was not supplied")
  if (missing(nSet)) stop("\'nSet\' was not supplied")
  if (missing(nRun)) stop("\'nRun\' was not supplied")

  require(reticulate)
  require(data.table)
  require(dplyr)

  reticulate::use_python(PyPath)
  # 22.08.26 added - "inst/"
  script <- paste(system.file(package="CINCIA"), "Detecting_Scr.py", sep = "/")
  source_python(script)

  cat("Detection result files will be saved in ", OutPath, "\n\n")
  OutPath = RunningScrublet(InPath, OutPath, type, nSet, nRun, estRate=est.rate)
}



#' Doublet scoring through Hybrid of scds, a R-based doublet detection program (Bais and Kostka, 2019)
#'
#' Doublet score (probability) of each cell in the input data or simulated data is yielded by Hybrid of Scds.
#'
#' @param InPath Path to the input directory in which input count matrix (for type='raw') or folders named 'simul_#' which contain
#' count matrices of doublet-simulated datasets (for type='simul') exist.
#' @param OutPath Path to the output directory in which the results files are saved.
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param nSet The number of dataset to run.
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @export
DoubletScoring_Hybrid <- function( InPath, OutPath, type, nSet, nRun=10){

  if (missing(InPath)) stop("\'InPath\' was not supplied")
  if (grepl("/$", InPath) == FALSE) InPath <- paste0(InPath, '/')
  if (missing(OutPath)) stop("\'OutPath\' was not supplied")
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  if (missing(type)) stop("\'type\' was not supplied")
  if (missing(nSet)) stop("\'nSet\' was not supplied")
  if (missing(nRun)) stop("\'nRun\' was not supplied")

  require(Seurat)
  require(scds)
  require(SingleCellExperiment)
  require(data.table)
  require(dplyr)

  #outScores.All <- c()
  start_time <- Sys.time()

  for (i in 1:nSet) {
    cat("===========================================================================================\n")
    cat("                                Dataset number ", i, "\n")
    cat('   Detection algorithm : Hybrid', "\n")

    # open expression matrix
    if (type=='raw') {
          cat("   Data type : Raw count matrix\n")
          setname <- paste('sample_', i, sep='')
          mtxFile <- list.files(InPath)[which(endsWith(list.files(InPath), "ExpMtx.nonZero.txt")==TRUE)]
          cat('   - Input file : ', paste0(InPath, mtxFile), "\n")

          #outdir <- OutPath
          outdir <- paste0(OutPath, '1.Raw/')
          if (dir.exists(outdir)==FALSE) dir.create(outdir)
          if (grepl("/$", outdir) == FALSE) outdir <- paste0(outdir, '/')

          cat('   - Output directory : ', outdir, "\n")
          cat("===========================================================================================\n")

          # 22.08.26 modified - fread
          Mtx <- fread(paste0(InPath, mtxFile), check.names = FALSE) %>% as.data.frame()
          rownames(Mtx) <- Mtx$V1
          Mtx <- Mtx %>% select(!V1)
          Mtx <- as.matrix(Mtx)

    } else {
          cat("   Data type : Simulated count matrix\n")
          setname <- paste0('simul_', i)
          indir <- paste0(InPath, setname, '/')
          mtxFile <- list.files(indir)[which(endsWith(list.files(indir), "simulDbls.ExpMtx.txt")==TRUE)]
          cat('   - Input file : ', paste0(indir, mtxFile), "\n")

          outdir <- paste0(OutPath, setname)
          if (dir.exists(outdir)==FALSE) dir.create(outdir)
          if (grepl("/$", outdir) == FALSE) outdir <- paste0(outdir, '/')
          cat('   - Output directory : ', outdir, "\n")
          cat("===========================================================================================\n")

          # 22.08.26 modified - fread
          Mtx <- fread(paste0(indir, mtxFile), check.names = FALSE) %>% as.data.frame()
          rownames(Mtx) <- Mtx$V1
          Mtx <- Mtx %>% select(!V1)
          Mtx <- as.matrix(Mtx)

          #Mtx <- read.delim(paste0(indir, mtxFile), header = T, sep='\t',  check.names = FALSE, row.names = 1)
          #Mtx <- as.matrix(Mtx)
    }

    sceobj <- SingleCellExperiment(list(counts=Mtx))
      # librarySizeFactor filtering
    sceobj <- sceobj[, colSums(counts(sceobj)) > 0]
    outScores <- c()

    for (j in 1:nRun){
        cat("=====================> Iteration ", j," <=====================\n")
        res <- cxds_bcds_hybrid(sceobj, verb = TRUE, force = TRUE)
        iter <- paste('iter_', j, sep='')
        outScores[[iter]] <- colData(res)
    }
    end_time <- Sys.time()
    #outScores.All[[setname]] <- outScores

    if (type=='raw') cat("\n--------------------> Finished <--------------------\n")
    else cat("\n--------------------> Dataset ", i, " Finished ! <--------------------\n\n")

    saveRDS(outScores, paste0( outdir, "Hybrid.res.Rds" ) )
  }
  cat("\nThe total execution time : ", TimeCheck(start_time, end_time) , "\n" )

}



#' Running scDblFinder (R-based doublet detection program)
#'
#' Doublet score (probability) of each cell in the input data or simulated data is yielded by scDblFinder.
#'
#' @param InPath Path to the input directory in which input count matrix (for type='raw') or folders named 'simul_#' which contain
#' count matrices of doublet-simulated datasets (for type='simul') exist.
#' @param OutPath Path to the output directory in which the results files are saved.
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param nSet The number of dataset to run.
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param est.rate Estimated doublet rates in your dataset.
#' @export
DoubletScoring_scDblFinder <- function( InPath, OutPath, type, nSet, nRun=10, est.rate=NULL){
  if (missing(InPath)) stop("\'InPath\' was not supplied")
  if (grepl("/$", InPath) == FALSE) InPath <- paste0(InPath, '/')
  if (missing(OutPath)) stop("\'OutPath\' was not supplied")
  if (grepl("/$", OutPath) == FALSE) OutPath <- paste0(OutPath, '/')
  if (missing(type)) stop("\'type\' was not supplied")
  if (missing(nSet)) stop("\'nSet\' was not supplied")
  if (missing(nRun)) stop("\'nRun\' was not supplied")

  require(scDblFinder)
  require(SingleCellExperiment)
  require(dplyr)
  require(data.table)

  #outScores.All <- c()
  start_time <- Sys.time()

  for (i in 1:nSet) {
    cat("===========================================================================================\n")
    cat("                                Dataset number ", i, "\n")
    cat('   Detection algorithm : scDblFinder', "\n")

    # open expression matrix
    if (type=='raw') {
        cat("   Data type : Raw count matrix\n")
        setname <- paste('sample_', i, sep='')
        mtxFile <- list.files(InPath)[which(endsWith(list.files(InPath), "ExpMtx.nonZero.txt")==TRUE)]
        cat('   - Input file : ', paste0(InPath, mtxFile), "\n")
        # outdir <- OutPath
        outdir <- paste0(OutPath, '1.Raw/')
        if (dir.exists(outdir)==FALSE) dir.create(outdir)
        if (grepl("/$", outdir) == FALSE) outdir <- paste0(outdir, '/')

        cat('   - Output directory : ', outdir, "\n")
        cat("===========================================================================================\n")

        # 22.08.26 modified - fread
        Mtx <- fread(paste0(InPath, mtxFile), check.names = FALSE) %>% as.data.frame()
        rownames(Mtx) <- Mtx$V1
        Mtx <- Mtx %>% select(!V1)
        Mtx <- as.matrix(Mtx)

        #Mtx <- read.delim(paste0(InPath, mtxFile), header = T, sep='\t',  check.names = FALSE, row.names = 1)
        #Mtx <- as.matrix(Mtx)
    } else {
        cat("   Data type : Simulated count matrix\n")
        setname <- paste0('simul_', i)
        indir <- paste0(InPath, setname, '/')
        mtxFile <- list.files(indir)[which(endsWith(list.files(indir), "simulDbls.ExpMtx.txt")==TRUE)]
        cat('   - Input file : ', paste0(indir, mtxFile), "\n")

        outdir <- paste0(OutPath, setname)
        if (dir.exists(outdir)==FALSE) dir.create(outdir)
        if (grepl("/$", outdir) == FALSE) outdir <- paste0(outdir, '/')
        cat('   - Output directory : ', outdir, "\n")
        cat("===========================================================================================\n")

        # 22.08.26 modified - fread
        Mtx <- fread(paste0(indir, mtxFile), check.names = FALSE) %>% as.data.frame()
        rownames(Mtx) <- Mtx$V1
        Mtx <- Mtx %>% select(!V1)
        Mtx <- as.matrix(Mtx)

        #Mtx <- read.delim(paste0(indir, mtxFile), header = T, sep='\t',  check.names = FALSE, row.names = 1)
        #Mtx <- as.matrix(Mtx)
    }

    sceobj <- SingleCellExperiment(list(counts=Mtx))
    sceobj <- sceobj[, colSums(counts(sceobj)) > 0]
    outScores <- c()

    for (j in 1:nRun){
      cat("=====================> Iteration ", j," <=====================\n")
      if ( is.null(est.rate)==FALSE ) {
        res <- scDblFinder::scDblFinder(sceobj, dbr = est.rate)
      } else{
        res <- scDblFinder::scDblFinder(sceobj)
      }
      iter <- paste('iter_', j, sep='')
      outScores[[iter]] <- colData(res)
    }
    end_time <- Sys.time()
    #outScores.All[[setname]] <- outScores

    if (type=='raw') cat("\n--------------------> Finished <--------------------\n")
    else cat("\n--------------------> Dataset ", i, " Finished ! <--------------------\n\n")

    saveRDS(outScores, paste0( outdir, "scDblFinder.res.Rds" ) )
  }
  cat("\nThe total execution time : ", TimeCheck(start_time, end_time) , "\n" )
}






#' Executing Doublet detection
#'
#' User-selected doublet detection algorithm among three, 'Scrublet', 'scDblFinder' and 'Hybrid', is executed to the user-selected dataset
#' which is either raw count matrix or simulated count matrices.
#'
#' @param PATH Path to the directory where raw/simulated count matrix and label-table exist.
#' Under this directory, a folder named '2.Simul/' already exists and '1.Raw/' will be generated when beginning doublet detection of raw data.
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param tool Doublet detection algorithm to execute. Select one among three algorithms( 'Scr'(Scrublet), 'scDF'(scDblFinder), 'Hyb'(Hybrid)).
#' Make sure you run all these tools for 'raw' and 'simul' data.
#' @param nSet The number of the dataset to run. Default to 1 for type 'raw' and 100 for type 'simul'.
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param PyPath Path to the Python executable in your computing environment (Python environment path).
#' Make sure the package 'Scrublet' is already installed.
#' @param est.rate Estimated doublet rates in your dataset, used in scDblFinder and Scrublet. If NULL, default value of each tool will be used.
#' @export
DoubletScoring <- function( PATH, type, tool, nSet, est.rate=NULL, nRun=10, PyPath=NULL){

  if (missing(PATH)) stop("\'PATH\' was not supplied")
  if (grepl("/$", PATH) == FALSE) PATH <- paste0(PATH, '/')

  if (missing(est.rate)) warning("\'est.rate\' was not supplied, using default of each detection tool.")
  if (missing(type)) stop('\'type\' was not supplied. Either \'raw\' or \'simul\' is required.')
  if (missing(tool)) stop("\'tool\' was not supplied.")
  if (missing(nSet)) {
    if (type=='raw'){
        warning('\'nSet\' was not supplied, using default 1 ')
        nSet = 1
    } else {
        warning('\'nSet\' was not supplied, using default 100')
        nSet = 100
    }
  }
  if (missing(nRun)) warning('\'nRun\' was not supplied, using default 10')
  if (missing(PyPath)){
      if (tool=='Scr') stop("\'PyPath\' was not supplied.")
  }

  if (type=="simul"){
    InPath <- paste0(PATH, "2.Simul/Input/")
    OutPath <- paste0(PATH, "2.Simul/Output/")
  } else{ #raw
    InPath <- PATH
    OutPath <- PATH
  }

  # Running each tool using data
  if (tool=='Hyb'){
      if (type=='raw') DoubletScoring_Hybrid( InPath=InPath, OutPath=OutPath, type='raw', nSet=nSet, nRun=nRun)
      else DoubletScoring_Hybrid( InPath=InPath, OutPath=OutPath, type='simul', nSet=nSet, nRun=nRun)
  } else {
      if (tool=='Scr'){
          if (type=='raw') DoubletScoring_Scrublet( PyPath, InPath, OutPath, type='raw', nSet=nSet, nRun=nRun, est.rate=est.rate )
          else DoubletScoring_Scrublet( PyPath, InPath, OutPath, type='simul', nSet=nSet, nRun=nRun, est.rate=est.rate)
      } else{
          if (type=='raw') DoubletScoring_scDblFinder( InPath=InPath, OutPath=OutPath, type='raw', nSet=nSet, nRun=nRun, est.rate=est.rate)
          else DoubletScoring_scDblFinder( InPath=InPath, OutPath=OutPath, type='simul', nSet=nSet, nRun=nRun, est.rate=est.rate)
      }
  }
}





#' Check Execution time
#'
#' @param start Execution start time
#' @param end Execution end time
#' @export
TimeCheck <- function(start, end) {
  s <- as.numeric(difftime(end, start, unit = "secs"))
  h <- floor(s / 3600)
  min <- floor((s - 3600 * h) / 60)
  sec <- s - 3600*h - 60*min
  paste0( sapply(c(h, min, sec),
                 function(x) { formatC(x, width = 2, format = "d", flag = "0")}), collapse = ":")
}
