% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwr_predict_par.R
\name{gwr_predict_par}
\alias{gwr_predict_par}
\title{GWR prediction}
\usage{
gwr_predict_par(fittedGWR = NULL, newdata, newcoords, longlat = NULL,
  se.fit = FALSE, ncores = NULL)
}
\arguments{
\item{fittedGWR}{A GWR object fitted with gwr_par()}

\item{newdata}{A new dataset used to predict the value (either data.frame of Spatial*Dataframe object)}

\item{newcoords}{A two-columns matrix with the coordinates of the new dataset if newdata is a dataframe Ignored otherwise}

\item{longlat}{TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers}

\item{se.fit}{Logical. TRUE if standard errors of the fit should be assessed.}

\item{ncores}{Number of cores in which the computation should be spanned}

\item{bandwith}{Bandwidth. Can be computed with gwr.sel.par()}
}
\value{
A SDF containing predicted values
}
\description{
This function predicts values of Y in location were it was not observed
}
