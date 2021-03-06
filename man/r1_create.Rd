% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/front_end.R
\name{r1_create}
\alias{r1_create}
\title{Create R1 knots on a grid a specificed number of knots wide}
\usage{
r1_create(locations, spatial_dimension, n_x = 10, buffer = 0)
}
\arguments{
\item{locations}{1d vector or 2d matrix of locations of the data.}

\item{spatial_dimension}{1 or 2}

\item{n_x}{how wide the grid should be in number of knots}

\item{buffer}{how much wider should the grid be than then data range}
}
\value{
vector or n by 2 matrix of knot locations
}
\description{
Create R1 knots on a grid a specificed number of knots wide
}
