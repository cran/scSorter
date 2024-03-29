% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities_temp.R
\name{xfindvariable_genes}
\alias{xfindvariable_genes}
\title{Select Highly Variable Genes}
\usage{
xfindvariable_genes(expr, ngenes = 2000)
}
\arguments{
\item{expr}{A matrix of input scRNA-seq data. Rows correspond to genes and columns correpond to cells.}

\item{ngenes}{The number of most variable genes to be selected.}
}
\value{
A vector of top highly variable genes with the total number determined by @ngenes option.
}
\description{
Select Highly Variable Genes following the vst approach. Please only use this function when you do not have access to Seurat package. More details are available in the vignette of this package.
}
