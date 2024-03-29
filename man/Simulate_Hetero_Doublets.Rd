% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Detect_Doublets.R
\name{Simulate_Hetero_Doublets}
\alias{Simulate_Hetero_Doublets}
\title{Heterotypic doublet simulation using reference set (singlets)}
\usage{
Simulate_Hetero_Doublets(refObj, nSimul = 100, nDbl, col.ref, OutPath)
}
\arguments{
\item{refObj}{Reference Set. A fully-processed Seurat object, containing putative singlets only.}

\item{nSimul}{The number of simulation runs. The number of simulated expression matrices as many as nSimul are created. Default to 100}

\item{nDbl}{The number of artificial doublets to generate for each simulation, which is recommended to be determined based on the number of cells
and multiplet rates of scRNA-seq method (ex. 10X Genomics Chromium) used. For each simulated dataset, artificial doublets within plus or minus 50 doublets
from nDbl are generated.}

\item{col.ref}{A Character string, name of the column which contains cell type label}

\item{OutPath}{Path to the output directory, where its subdirectories contain log files generated in previous steps.
New folders named 'simul_#'(total nSimul folders) are created and each simulated dataset is saved in there.}
}
\description{
Artificial heterotypic doublets are simulated using cells of reference Set, and combined with singlets together
into a gene expression matrix ( genes X cell barcodes ) and label table containing both singlets and doublets is generated.
}
