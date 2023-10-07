\name{tmda}


\description{
  f the corresponding normal distribution equal to \code{sd_norm}.
}

\name{tmda}
\alias{tmda}
\title{The TMDA Model}
\description{
  Fit a TMDA model.
}
\usage{
tmda(x_array, dims, y, subnum, cov.str = "general", ...)
}
\arguments{
  \item{x_array}{A matrix where each row is a vectorized (in column major order) tensor predictor.}
  \item{dims}{A vector contains the dimension of a tensor predictor.}
  \item{y}{A vector of class labels which has to be coded as continuous positive integers starting from 1.}
  \item{subnum}{A vector contains the number of mixture components for all classes.}
  \item{cov.str}{The covariance structure of fitted model. Possible options are "diag", "general", "scalar". Default is "general".}
}
\value{
  A list with the following components:
  \item{mu_list}{A list of centers.}
  \item{sig}{A list of square root of inverse covariance matrices (not covariance matrices).}
  \item{pi_list}{A list of prior probability of each mixture component.}
  \item{gamma}{A list of posterior probability of each mixture component.}
  \item{loglike}{The \eqn{log}-likelihood of the fitted model.}
}
\examples{
data(gsa)
x = gsa$x2
y = gsa$y2
cv_result = cv.tmda(x, c(8,16), y, 2:10, k=10)
min(cv_result$cv_error)
best_subnum = cv_result$`best subnum`
tmda_mod = tmda(x, c(8,16), y, rep(best_subnum,6))
tmda_result = predict.EM(x, tmda_mod)
tmda_pred = tmda_result$pred
tmda_err = loss(y, tmda_pred2)
tmda_err
}
\seealso{
  \code{\link{predict.EM}}, \code{\link{cv.tmda}}
}
\keyword{TMDA}
\keyword{classification}
\keyword{mixture models}
\keyword{EM algorithm}