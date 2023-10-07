tucker_prod = function(tnsr, mat_list, skip = NULL){
  dims = dim(tnsr)
  nlen = prod(dims)
  num_mat = length(mat_list)
  mat_param = makevec(mat_list)
  if(length(dims) != num_mat) stop("number of matrices not match tensor rank")
  if(any(dims != mat_param[[2]][2,])) stop("dimension not match")
  new_dims = mat_param[[2]][1,]
  output = array(0.0, new_dims)
  newlen = prod(mat_param[[2]][1,])
  if(is.null(skip) == TRUE){
    skip = rep(-1, length(dims))
    output = .Fortran("tucker_prod_f", as.numeric(output), 
                      as.integer(new_dims),
                      as.numeric(tnsr), as.integer(dims), 
                      as.integer(length(dims)), as.numeric(mat_param[[1]]), 
                      as.integer(length(mat_param[[1]])), skip)[[1]]
  }else{
    output = .Fortran("tucker_prod_f", as.numeric(output), 
                      as.integer(new_dims),
                      as.numeric(tnsr), as.integer(dims), 
                      as.integer(length(dims)), as.numeric(mat_param[[1]]), 
                      as.integer(length(mat_param[[1]])),
                      as.integer(skip), skip)[[1]]
  }
  dim(output) = new_dims
  return(output)
}

makevec = function(mat_list){
  output = list()
  output[[1]] = unlist(lapply(mat_list, as.vector))
  output[[2]] = sapply(mat_list, dim)
  return(output)
}


loss = function(y, pred){
  sum(y != pred) / length(y)
}