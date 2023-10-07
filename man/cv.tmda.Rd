\name{cv.tmda}
\alias{cv.tmda}
\title{Perform k-fold Cross-Validation on the TMDA Model}
\description{
  Performs k-fold cross-validation on a tensor mixture discriminant analysis model. It also iteratively fits the model with different numbers of mixture components to find the optimal number. This function assumes the number of mixture components is the same for each class.
}
\usage{
cv.tmda(x, dims, y, subnum = c(2, 3), cov_str = "general", k = 2)
}
\arguments{
  \item{x}{A matrix where each row is a vectorized (in column major order) tensor predictor.}
  \item{dims}{A vector containing the dimensions of tensor predictors.}
  \item{y}{A vector of class labels which must be coded as continuous positive integers starting from 1.}
  \item{subnum}{A vector or a list. When a vector is used, a vector of the numbers of mixture components is considered (assume the same number for each class). Default is c(2, 3). When a list is used, the number of mixture component can be different for each class.}
  \item{cov_str}{The covariance structure of the fitted model. Possible options are 'diag', 'general', and 'scalar'. Default is 'general'.}
  \item{k}{The number of folds in the k-fold cross-validation. Default is 2.}
}
\value{
  \item{cv_error}{A vector of the mean cross-validation error for each number of mixture components.}
  \item{best subnum}{The number of mixtures with the smallest cross-validation error.}
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
\code{\link{tmda}}, \code{\link{predict.EM}}
}

\keyword{cross-validation}
\keyword{classification}
\keyword{mixture models}
\keyword{EM algorithm}
