\name{csa}
\alias{csa}
\docType{data}
\title{
Colorimetric sensor array (CSA) data}
\description{
The dataset comprises observations collected from CSA experiments involving volatile chemical toxicants (VCT). Chemical dyes were exposed to VCT and the resulting colors of the dyes were recorded to determine the class of the VCT.

}
\usage{
data(csa)
}

\format{
  A data matrix \emph{x} contains 294 colorimetric sensor array observations. Each observation represents a \eqn{36\times3} dimensional array being vectorized in a column-major format. A numeric vector \emph{y} contains labels. 
  \describe{
    \item{\code{x}}{Colorimetric sensor array matrix.}
    \item{\code{y}}{Class label ranging from 1 to 21.}
  }
}


\details{
This dataset is provided in the Supplementary matrial of Zhong (2015). The data consists of two concentrations. We reorganize the colorimetric sensor arrays into a matrix. Originally, each row of a colorimetric sensor array represents a dye in red, green and blue.
}
\source{
Wenxuan Zhong and Kenneth S. Suslick (2015). "Matrix discriminant analysis with application to colorimetric sensor array data" \emph{Technometrics} \bold{57}(4), 524--534.
}
\keyword{datasets}
