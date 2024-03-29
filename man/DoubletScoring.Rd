% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Detect_Doublets.R
\name{DoubletScoring}
\alias{DoubletScoring}
\title{Executing Doublet detection}
\usage{
DoubletScoring(
  PATH,
  type,
  tool,
  nSet,
  est.rate = NULL,
  nRun = 10,
  PyPath = NULL
)
}
\arguments{
\item{PATH}{Path to the directory where raw/simulated count matrix and label-table exist.
Under this directory, a folder named '2.Simul/' already exists and '1.Raw/' will be generated when beginning doublet detection of raw data.}

\item{type}{Type of the input dataset, either 'raw' or 'simul'.}

\item{tool}{Doublet detection algorithm to execute. Select one among three algorithms( 'Scr'(Scrublet), 'scDF'(scDblFinder), 'Hyb'(Hybrid)).
Make sure you run all these tools for 'raw' and 'simul' data.}

\item{nSet}{The number of the dataset to run. Default to 1 for type 'raw' and 100 for type 'simul'.}

\item{est.rate}{Estimated doublet rates in your dataset, used in scDblFinder and Scrublet. If NULL, default value of each tool will be used.}

\item{nRun}{An integer indicating the number of running iteration with different random seeds for each dataset. Default to 10}

\item{PyPath}{Path to the Python executable in your computing environment (Python environment path).
Make sure the package 'Scrublet' is already installed.}
}
\description{
User-selected doublet detection algorithm among three, 'Scrublet', 'scDblFinder' and 'Hybrid', is executed to the user-selected dataset
which is either raw count matrix or simulated count matrices.
}
