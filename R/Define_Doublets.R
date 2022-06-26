
#' Binding the doublet detection results
#'
#' Doublet scores calculated from three tools (Hybrid, scDblFinder, Scrublet ) are loaded.
#' For each tool and each dataset (either raw or simulated), detection results (.txt) are read and bound,
#' resulting in a list of data frame containing three doublet scores of each cell.
#'
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param nSet The number of the dataset to run. Default to 1 for type 'raw' and 100 for type 'simul'.
#' @param InPath Path to the input directory in which input count matrix (for type='raw') or folders named 'simul_#' which contain
#' count matrices of doublet-simulated datasets (for type='simul') exist.
#' @return A list of doublet scores calculated from three methods for each cell in the data.
#' @export
Bind_Doublet_Scores <- function( type, nRun=10, nSet, InPath ) {
  if (missing(type))
      stop('\'type\' was not supplied. either \'raw\' or \'simul\' for \'type\' is required.')
  if (missing(nRun))
      warning('\'nRun\' was not supplied, using default 10')

  if (missing(nSet)) {
      if (type=='raw'){
          warning('\'nSet\' was not supplied, using default 1 ')
          nSet = 1
      } else {
          warning('\'nSet\' was not supplied, using default 100')
          nSet = 100
      }
  }
  if (missing(InPath))
      stop('\'InPath\' was not supplied, Path to the Input file directory is required.')

  if (grepl("/$", InPath) == FALSE) { InPath <- paste0(InPath, '/') }

  DblScores.All <- list()
  Hybrid.scores <- list()
  Scrublet.scores <- list()
  scDblFinder.scores <- list()

  # 1. Hybrid (scds)

  cat("==== 1) Reading Hybrid results ====\n\n")
  for (i in 1:nSet) {
    if (type=="raw") {
        cat("  ## raw dataset ", i, '\n')
        set_name <- paste('data_', i, sep='')
    } else {
        cat("  ## simulation set ", i, '\n')
        set_name <- paste('simul_', i, sep='')
    }

    # Read in the doublet scores
    for(j in 1:nRun){
      if (type=="raw") {
          file <- paste0(InPath, "1.Raw/Hybrid.res.Rds")
          if (file.exists(file)==FALSE){
              file <- paste0(InPath, "1.Raw/ Hybrid.res.Rds")
          }
      } else {
          file <- paste0(InPath, "2.Simul/Output/", set_name, "/Hybrid.res.Rds")
          if (file.exists(file)==FALSE){
              file <- paste0(InPath, "2.Simul/Output/", set_name, "/ Hybrid.res.Rds")
          }
      }
      result <- readRDS(file)

      if(j==1){ res <- data.frame( row.names = rownames(result[[j]]), result[[j]]$hybrid_score) }
      else{     res <- cbind(res, result[[j]]$hybrid_score)} }

    res <- cbind(res, as.vector(apply(res, MARGIN=1, mean)))
    colnames(res) <- c(rep(paste('Hyb.iter', 1:nRun, sep='')), 'Hyb.Average')
    Hybrid.scores[[i]] <- as.data.frame(res)
  }


  # 2. scDblFinder
  cat("\n==== 2) Reading scDblFinder results ====\n\n")
  for(i in 1:nSet){
    if (type=="raw") {
      cat("  ## raw dataset ", i, '\n')
      set_name <- paste('data_', i, sep='')
    } else {
      cat("  ## simulation set ", i, '\n')
      set_name <- paste('simul_', i, sep='')
    }

    # Read in the doublet scores
    for(j in 1:nRun){
      if (type=="raw") {
          file <- paste0(InPath, "1.Raw/scDblFinder.res.Rds")
      } else {
          file <- paste0(InPath, "2.Simul/Output/",  set_name, "/scDblFinder.res.Rds")
      }
      result <- readRDS(file)

      if(j==1){ res <- data.frame( row.names = rownames(result[[j]]), result[[j]]$scDblFinder.score) }
      else{     res <- cbind(res, result[[j]]$scDblFinder.score)} }

    res <- cbind(res, as.vector(apply(res, MARGIN=1, mean)))
    colnames(res) <- c(rep(paste('scDF.iter', 1:nRun, sep='')), 'scDF.Average')
    scDblFinder.scores[[i]] <- as.data.frame(res)
  }


  # 3. Scrublet

  cat("\n==== 3) Reading Scrublet results ====\n\n")
  for (i in 1:nSet) {
    if (type=="raw") {
        cat("  ## raw dataset ", i, '\n')
        set_name <- paste0('data_', i)
    } else {
        cat("  ## simulation set ", i, '\n')
        set_name <- paste0('simul_', i)
    }

    # collect doublet scores
    for ( j in 1:nRun ){
      if (type=="raw") {
          file <- paste0( "1.Raw/scr_", j,".scrublet_result.txt")
      } else {
          file <- paste(  "2.Simul/Output/", set_name, "/scr_", j,".scrublet_result.txt", sep='')
      }
      result <- read.delim(paste0(InPath, file), header=T, sep='\t')

      if (j==1) {
          res <- data.frame( row.names = result$Cell_barcode, result$Scrublet_score)
      } else {
          res <- cbind(res, result$Scrublet_score)
      }
    }

    res <- cbind(res, as.vector(apply(res, MARGIN=1, mean)))
    colnames(res) <- c(rep(paste('Scr.iter', 1:nRun, sep='')), 'Scr.Average')
    Scrublet.scores[[i]] <- as.data.frame(res)
  }

  cat("\n Finished  \n\n")

  if (type=="raw") name = 'data_'
  else name = 'simul_'

  names(Hybrid.scores) <- rep(paste0(name, 1:nSet))
  names(scDblFinder.scores) <- rep(paste0(name, 1:nSet))
  names(Scrublet.scores) <- rep(paste0(name, 1:nSet))

  DblScores.All <- list( "Hybrid"= Hybrid.scores,
                         "scDblFinder" = scDblFinder.scores,
                         "Scrublet" = Scrublet.scores)
  return(DblScores.All)
}





#' Normalizing doublet scores between 0 and 1
#'
#' For each dataset and detection tool, doublet scores previously calculated in 'Bind_Doublet_Scores()'
#' are averaged and normalized to range (0,1) by min-max scaling.
#'
#' @param type Type of the input dataset, either 'raw' or 'simul'.
#' @param dblobj A list of doublet scores, resulted from function Bind_Doublet_Scores().
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param nSet The number of datasets to run. Default to 1 for 'raw' and 100 for 'simul'.
#' @export
Scale_AvgScores <- function( type, dblobj, nRun=10, nSet){

  if (missing(type))
      stop('\'type\' was not supplied. Type of the dataset, either \'raw\' or \'simul\' is required.')
  if (missing(dblobj))
     stop('\'dblobj\', object resulted from function \'Bind_Doublet_Scores\' is required.')
  if (missing(nRun)){
     warning('\'nRun\' was not supplied, using default 10')
  }
  if (missing(nSet)){
      if (type=='raw'){
          warning('\'nSet\' was not supplied, using default 1 ')
          nSet = 1
      } else {
          warning('\'nSet\' was not supplied, using default 100')
          nSet = 100
      }
  }

  result <- c()
  cat('   # min-max scaling of average scores\n')

  for (i in 1:length(names(dblobj))) {
    tool <- names(dblobj)[i]
    cat('     -', i, tool, '\n')

    for (j in 1:nSet){
        setname <- names(dblobj[[tool]])[j]
        avg <- dblobj[[tool]][[j]][nRun+1]
        avg.scale <- as.data.frame( apply(avg, MARGIN = 2, FUN=function(x){ (x-min(x))/diff(range(x)) }) )

        if (i==1){
            # check rownames
            rownames(avg.scale) <- gsub("\\.", "-", rownames(avg.scale))
            result[[setname]] <- avg.scale
        } else {
            result[[setname]] <- cbind(result[[setname]],
                                   avg.scale[ rownames(result[[setname]]), 1 ] )
            if (i==length(names(dblobj))){
                names(result[[setname]]) <- names(dblobj)
                # Assign class label to each cell
                if (type=='simul'){
                    result[[setname]]['Actual.type'] <- rep('singlet', nrow(result[[setname]]))
                    result[[setname]]['Actual.type'] <- ifelse( startsWith(rownames(result[[setname]]), 'dbl_'),
                                                                'doublet', 'singlet' )
          }
        }
      }
    }
  }
  return(result)
}




#' Defining doublets using normalized doublet scores
#'
#' For each doublet detection method, probability density functions (PDFs) of singlets and doublets in the simulated data
#' are calculated and then putative doublets in the input data are defined using Gaussian naive bayes classification, based on those PDFs.
#'
#' @param SimulSet A list containing normalized doublet scores of simulated datasets, resulted from function 'Scale_AvgScores'.
#' @param RawSet A data-frame containing normalized doublet scores of input (raw) dataset, resulted from function 'Scale_AvgScores'.
#' @export
Doublet_calling <- function( SimulSet, RawSet ){

  if (missing(SimulSet))
    stop('\'SimulSet\' was not supplied. ')
  if (missing(RawSet))
    stop('\'RawSet\' was not supplied.')

  require(dplyr)
  require(reshape2)

  data <- bind_rows(SimulSet)

  # 1. Calculate mean and sd according to each class (singlet/doublet) and each tool

  means <- c(); sds <- c(); tools <- c(); classes <- c()
  all.tools <- setdiff(names(SimulSet[[1]]), c("Actual.type"))

  for (tool in all.tools){
      for (class in c('singlet','doublet')){
          means <- c(means, mean(data[data$Actual.type==class, tool]))
          sds <- c(sds, sd(data[data$Actual.type==class, tool]))
          tools <- c(tools, tool)
          classes <- c(classes, class)
      }
  }
  stats <- melt( data.frame( 'mean' = means, 'sd'=sds, 'tool'=tools, 'class'=classes ),
                 id.vars = c('class', 'tool'))
  print(stats)
  cat("\n\n")


  # 2. Calculate Prior probability based on simulation data

  ndbl <- c(); nsgl <- c()
  for (i in 1:length(SimulSet)) {
      ndbl <- c(ndbl, table(SimulSet[[i]]['Actual.type'])[['doublet']])
      nsgl <- c(nsgl, table(SimulSet[[i]]['Actual.type'])[['singlet']])
  }
  prior.Dbl <- mean(ndbl) / sum( mean(nsgl), mean(ndbl) )
  prior.Sgl <- mean(nsgl) / sum( mean(nsgl), mean(ndbl) )

  cat(" Estimated number of doublets = ", mean(ndbl), "\n")
  cat(" Estimated number of singlets = ", mean(nsgl), "\n")
  cat(" ===> Prior probability of doublets and singlets : ", prior.Dbl, ", ", prior.Sgl , "\n\n")


  # 3. posterior probabilities of each cell in test set (raw data)

  log.Posterior <- data.frame( row.names <- c(), post.Doublet <- c(), post.Singlet <- c() )

    # for each cell..
  Rawdf <- RawSet[[1]]
  for (i in 1:nrow(Rawdf)){
      cell <- rownames(Rawdf)[i]
      Doublet <- c(); Singlet <- c()

      # calculate posterior probabilities
      # for class
      for( class in c('singlet', 'doublet')){
        # for detection method
          for(j in 1:length(all.tools)){
              tool <- all.tools[j]
              mean <- stats[stats$class==class & stats$tool==tool & stats$variable=='mean', 'value']
              sd <- stats[stats$class==class & stats$tool==tool & stats$variable=='sd','value']

              if (class=="singlet") {
                  Singlet[j] <- log( dnorm( Rawdf[i,tool], mean=mean, sd=sd ))
              } else {
                  Doublet[j] <- log( dnorm( Rawdf[i,tool], mean=mean, sd=sd ))
              }
          }
      }

    # obtain log-transformed doublet score
    post.dbl <- log(prior.Dbl) + sum(Doublet)
    post.sgl <- log(prior.Sgl) + sum(Singlet)
    pred <- ifelse( post.dbl > post.sgl, 'doublet', 'singlet')
    log.Posterior <- rbind(log.Posterior,
                           data.frame(row.names=cell, post.Doublet=post.dbl, post.Singlet=post.sgl, predict=pred))
  }
  return(log.Posterior)
}




#' Executing all the functions for defining putative doublets in the data
#'
#' In this function, putative doublets exist in the data are defined.
#'
#' @param nRun An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10
#' @param nSet_s The number of simulated datasets. Default to 100.
#' @param nSet_a The number of actual input (raw) dataset. Default to 1
#' @param PATH Path to the input directory in which input count matrix (for type='raw') or folders named 'simul_#' which contain
#' count matrices of doublet-simulated datasets (for type='simul') exist.
#' @param PATH_out Path to the directory to save result of each step
#' @export
DefineDoublets <- function( nRun, nSet_s, nSet_a, PATH, PATH_out ){

  if (missing(nSet_s)){
      warning('\'nSet_s\' was not supplied, using default 100')
      nSet_s = 100
  }
  if (missing(nSet_a)){
      warning('\'nSet_a\' was not supplied, using default 1 ')
      nSet_a = 1
  }
  if (missing(nRun)){
      warning('\'nRun\' was not supplied, using default 10')
  }

  if (missing(PATH)) stop('\'PATH\' was not supplied')
  if (missing(PATH_out)) stop('\'PATH_out\' was not supplied')

  if (grepl("/$", PATH) == FALSE) PATH <- paste0(PATH, '/')
  if (grepl("/$", PATH_out) == FALSE) PATH_out <- paste0(PATH_out, '/')

  Outdir <- paste0(PATH_out, '2.Doublet/')
  dir.create(Outdir)

  # 1. Combining the doublets

  cat("\nCombining doublet detection results of simulated data.....\n\n")
  simul.DblResult <- Bind_Doublet_Scores( type='simul', nRun=10, nSet=nSet_s, PATH )

  cat("\nCombining doublet detection results of raw data.....\n\n")
  raw.DblResult <- Bind_Doublet_Scores( type='raw',  nRun=10, nSet=nSet_a,  PATH )

  saveRDS( list('simul'=simul.DblResult, 'raw'=raw.DblResult),
           paste0(Outdir, "1-1.Dblresult.All.Rds") )

  # 2. Scaling doublet scores

  cat("\n\nMin-max scaling of doublet scores.....\n")
  cat("\n 1) Simulated data\n")
  sc.simul <- Scale_AvgScores( type='simul', dblobj=simul.DblResult, nRun=10, nSet=nSet_s  )

  cat("\n 2) raw data\n")
  sc.raw <- Scale_AvgScores( type='raw', dblobj=raw.DblResult, nRun=10, nSet=nSet_a )
  saveRDS( list('Simul'=sc.simul, 'Raw'=sc.raw),
           paste0(Outdir, "1-2.Dblresult.All.scaled.Rds"))

  # 3. Gaussian Naive bayes classification

  cat("\n\nClassifying putative doublets.....\n")
  NB.result <- Doublet_calling(sc.simul, sc.raw)
  saveRDS(NB.result, paste0(Outdir, "2.NB.Putative.Doublets.Rds"))

  return(NB.result)
}
