\name{gsa}
\alias{gsa}
\docType{data}
\title{
Gas sensor array (GSA) data}
\description{
The dataset contains measurements from 16 chemical sensors exposed to 6 gases. Sensor responses are converted into real scalars resulting in a 128-dimensional feature vector (8 features x 16 sensors).

}
\usage{
data(gsa)
}

\format{
  Data matrices \emph{x1} and \emph{x2} contain 3613 and 470 gas sensor array observations respectively. Numeric vectors \emph{y1} and \emph{y2} contain labels correspondingly. Each sensor array represents a \eqn{8\times16} dimensional array being vectorized in a column-major format.  
  \describe{
    \item{\code{x1}}{Gas sensor array matrix.}
    \item{\code{x2}}{Gas sensor array matrix.}
    \item{\code{y1}}{Class label ranging from 1 to 6.}
    \item{\code{y2}}{Class label ranging from 1 to 6.}
  }
}


\details{
This dataset is provided in Vergara et al. (2012). The original data is organized into batches according to the month of collection. Two batches are given here as examples of large sample size (3613) and small sample size (470), and the concentration level information is omitted.
}
\source{
Alexander Vergara, Shankar Vembu, Tuba Ayhan, Margaret A. Ryan, Margie L. Homer, Ramón Huerta (2012). "Chemical gas sensor drift compensation using classifier ensembles" \emph{Sensors and Actuators B: Chemical} Volumes(166–167), 320-329.
}
\keyword{datasets}
