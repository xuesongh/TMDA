\name{tucker_prod}
\alias{tucker_prod}
\title{Tucker Product of a Tensor and a List of Matrices}
\description{
  This function computes the tucker product of a tensor and a list of matrices. It returns an array that is the result of the tucker product operation.
}
\usage{
tucker_prod(tnsr, mat_list, skip = NULL)
}
\arguments{
  \item{tnsr}{An array representing the tensor to be multiplied.}
  \item{mat_list}{A list of matrices to multiply with the tensor. The number of matrices should match the rank of the tensor, and the dimensions of the matrices should match the corresponding dimensions of the tensor.}
  \item{skip}{An optional argument. If provided, it should be an integer that indicates a dimension to be skipped in the tucker product operation.}
}
\value{
  Returns an array that is the result of the tucker product operation.
}

\examples{
# Create a tensor
tnsr <- array(1:8, c(2, 2, 2))

# Create a list of matrices
mat_list <- list(matrix(1:4, 2), matrix(1:4, 2), matrix(1:4, 2))

# Compute the tucker product
result <- tucker_prod(tnsr, mat_list)

print(result)
}
\keyword{tucker product}
