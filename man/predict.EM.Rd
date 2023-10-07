\name{predict.EM}
\title{Predict function for the TMDA model}
\description{
This function predicts class labels for new data based on the output from the \code{tmda} function.
}
\usage{
predict.EM(newx, output)
}
\arguments{
  \item{newx}{A matrix of test data where each row is a vectorized (in column major order) tensor predictor.}
  \item{output}{A list of \code{tmda} output.}
}
\value{
  \item{prob}{A matrix of probabilities where each row is the posterior probability of each class.}
  \item{pred}{A vector of predicted class labels.}
}
\seealso{
  \code{\link{tmda}}
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
\keyword{prediction}
\keyword{classification}
\keyword{mixture models}
\keyword{EM algorithm}
