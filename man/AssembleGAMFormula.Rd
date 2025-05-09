% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AssembleGAMFormula.R
\name{AssembleGAMFormula}
\alias{AssembleGAMFormula}
\title{Make a GAM formula}
\usage{
AssembleGAMFormula(yvar, gam.table, hgam = F)
}
\arguments{
\item{yvar}{Name of dependent variable for gam models}

\item{gam.table}{Data frame of parameters for GAM formula}

\item{hgam}{Logical; do you want an hgam formula}
}
\value{
Returns a formula object, or list of formulas for hgam
}
\description{
Improved version designed to use the tables produced by the Autodetect functions
}
