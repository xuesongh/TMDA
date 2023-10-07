predict.EM = function(newx, output){
  mu_list = output[["mu"]]
  pi_list = output[["pi"]]
  sig = output[["sig"]]
  PI = output[["PI"]]
  cov_str = output["cov_str"]
  dims = output[["dims"]] ###
  n = nrow(newx) ###
  nclass = length(PI)
  subnum = sapply(pi_list, length)
  prob = matrix(0, n, nclass)
  if(cov_str == "diag"){
    sig = unlist(lapply(sig, diag))
    output = .Fortran("pred_diag", as.numeric(t(newx)), 
                      as.integer(dims), as.integer(length(dims)),
                      as.numeric(prob), as.integer(n), as.integer(nclass),
                      as.integer(subnum), as.numeric(unlist(mu_list)), 
                      as.numeric(sig), as.numeric(unlist(pi_list)),
                      as.numeric(PI))
  }else{
    output = .Fortran("pred_f", as.numeric(t(newx)), 
                      as.integer(dims), as.integer(length(dims)),
                      as.numeric(prob), as.integer(n), as.integer(nclass),
                      as.integer(subnum), as.numeric(unlist(mu_list)), 
                      as.numeric(unlist(sig)), as.numeric(unlist(pi_list)),
                      as.numeric(PI))
  }

  prob_mat = matrix(output[[4]], n, nclass)
  pred = apply(prob_mat, 1, which.max)
  out = list()
  out[["prob"]] = (prob_mat)
  out[["pred"]] = pred
  return(out)
}

EM_fortran = function(x_array, dims, y, subnum, mu_list, pi_list = NULL, 
                      sig = NULL, gamma = NULL, cov_str = "diag",
                      maxit = 100, tol=1e-5){
  loglike = rep(-1, maxit)
  loglike[maxit] = tol
  lab = sort(unique(y))
  nclass = length(lab)
  if(length(subnum) != nclass) stop("numbers of subclasses are not correctly specified")
  if(length(mu_list) != nclass) stop("number of label and initial mu not match")
  # create gamma list
  glen = sum(subnum[y])
  if(is.null(gamma)){
    gamma = numeric(glen)
  }
  # initialize pi for each subclass by uniform
  if(is.null(pi_list)){
    pi_list = list()
    for (i in 1:nclass) {
      pi_list[[i]] = rep(1/subnum[i], subnum[i])
    }
  }
  # initialize sigma by identity matrix
  n = nrow(x_array)
  nlen = ncol(x_array)
  M = length(dims)

  if(cov_str == "diag"){
    if(is.null(sig)){
      sig = rep(1, sum(dims))
    }
    output = .Fortran("em_diag", as.numeric(t(x_array)), 
                as.integer(dims), as.integer(length(dims)),
                as.integer(y), as.integer(n), as.integer(nclass),
                as.integer(subnum), as.numeric(unlist(mu_list)), 
                as.numeric(unlist(sig)), as.numeric(unlist(pi_list)),
                as.numeric(unlist(gamma)), as.integer(glen),
                as.integer(maxit), as.numeric(loglike))
    # rearrange sig
    sig = list()
    index = 0
    for(m in 1:M){
      sig[[m]] = diag(output[[9]][(index+1) : (index+dims[m])])
      index = index + dims[m]
    }
  }else{
    if(is.null(sig)){
      sig = list()
      for (i in 1:M) {
        sig[[i]] = diag(dims[i])
      }
    }
    output = .Fortran("em", as.numeric(t(x_array)), 
                    as.integer(dims), as.integer(length(dims)),
                    as.integer(y), as.integer(n), as.integer(nclass),
                    as.integer(subnum), as.numeric(unlist(mu_list)), 
                    as.numeric(unlist(sig)), as.numeric(unlist(pi_list)),
                    as.numeric(unlist(gamma)), as.integer(glen),
                    as.integer(maxit), as.numeric(loglike))
    # rearrange sig
    index = 0
    for(m in 1:M){
      sig[[m]] = matrix(output[[9]][(index+1) : (index+dims[m]**2)], dims[m], dims[m])
      index = index + dims[m] ** 2
    }
  }
  PI = table(y) / sum(table(y))
  out = list()
  gamma_list = list()
  record = numeric(nclass)
  for (j in 1:nclass) {
    index = sum(subnum[1:j]) - subnum[j]
    pi_list[[j]] = output[[10]][(index+1):(index+subnum[j])]
    for(r in 1:subnum[j]){
      index = (sum(subnum[1:j]) - subnum[j] + r-1)*nlen
      mu_list[[j]][[r]] = array(output[[8]][(index+1):(index+nlen)], dims)
    }
  }
  index = 0
  for(i in 1:n){
    if(record[y[i]] == 0){
      gamma_list[[y[i]]] = output[[11]][(index+1):(index+subnum[y[i]])]
    }else{
      gamma_list[[y[i]]] = rbind(gamma_list[[y[i]]], 
                                 output[[11]][(index+1):(index+subnum[y[i]])])
    }
    record[y[i]] = record[y[i]] + 1
    index = index + subnum[y[i]]
  }
  out[["mu"]] = mu_list
  out[["pi"]] = pi_list
  out[["sig"]] = sig
  out[["gamma"]] =  gamma_list
  out[["PI"]] = as.vector(PI)
  out[["cov_str"]] = cov_str
  out[["loglike"]] = output[[14]][output[[13]]]
  out[["iteration"]] = output[[13]]
  out[["dims"]] = dims
  out
}

tmda = function(x_array, dims, y, subnum, cov.str = "general", maxit = 50, 
                tol = 1e-5, mu = NULL, pi_list = NULL,
                sig = NULL, gamma = NULL){
  ##initialize mu by k-means
  M = length(dims)
  nclass = length(table(y))
  mask = NULL; grand_mean = NULL
  if(length(dim(x_array)) != 2) stop("incorrect x dimension")
  if(dim(x_array)[2] != prod(dims)) stop("incorrect dims")
  if(length(subnum) == 1) subnum = rep(subnum, nclass)
  if(is.null(mu)){
    k_means_list = list()
    mu = list()
    pi_list = list()
    for (j in 1:length(subnum)) {
      k_means_list[[j]] = kmeans(x_array[y == j,], subnum[j])
      pi_list[[j]] = k_means_list[[j]]$size / sum(k_means_list[[j]]$size)
      mu[[j]] = list()
      for(r in 1:subnum[j]){
        mu[[j]][[r]] = array(k_means_list[[j]]$centers[r,], dims)
      }
    }
    #print("kmeans initial finished")
  }
  if(cov.str == "diag"){
    out = EM_fortran(x_array, dims, y, subnum, mu, pi_list, sig, gamma, 
                     cov_str = "diag", maxit, tol)
  }else{
    out = EM_fortran(x_array, dims, y, subnum, mu, pi_list, sig, gamma,
                     cov_str = "general", maxit, tol)
  }
  out[["mask_level"]] = qt
  return(out)
}

cv.tmda <- function(x, dims, y, subnum = c(2, 3), cov_str = "general", k = 2,
                    iters = 500){
  n <- length(y)
  foldid <- sample(rep(1:k, length = n))
  cv_error <- matrix(-1, nrow = length(subnum), ncol = k)
  countj <- 0
  iterations <- matrix(-1, nrow = length(subnum), ncol = k)
  inits <- list()
  mask = NULL
  for (j in subnum) {
    print(paste("CV number of mixtures: ", toString(j)))
    countj <- countj + 1
    for (i in 1:k) {
      test_index <- which(foldid == i)
      iteration <- 0
      if (i == 1) {
        dif <- -1
        count <- 1
        er <- c()
        er[1] <- -1
        while (dif < 0 | count < 2) {
          if (count == 1) {
            if (all(j == subnum[[1]])) {
              mod0 <- tmda(x[-test_index,], dims, y[-test_index],
                j, cov_str,
                maxit = ceiling(iters / 10),
                tol = 1e-10
              )
            } else {
              mod0 <- tmda(x[-test_index,], dims, y[-test_index],
                j, cov_str,
                maxit = ceiling(iters / 10),
                tol = 1e-10
              )
            }
          } else {
            mod0 <- tmda(x[-test_index,], dims, y[-test_index],
              j, cov_str,
              maxit = ceiling(iters / 10), tol = 1e-10,
              mu = mod0$mu, pi_list = mod0$pi,
              sig = mod0$sig, gamma = mod0$gamma
            )
          }
          count <- count + 1
          iteration <- iteration + ceiling(iters / 10)
          pred <- predict.EM(x[test_index,], mod0)$pred
          er[count] <- loss(y[test_index], pred)
          dif <- er[count] - er[count - 1]
        }
        inits[[countj]] <- list()
        inits[[countj]][["mu"]] <- mod0$mu
        inits[[countj]][["pi"]] <- mod0$pi
        inits[[countj]][["sig"]] <- mod0$sig
        inits[[countj]][["gamma"]] <- mod0$gamma
        err <- er[count]
        iters <- max(1, iters * (count - 2) / 10)
      } else {
        mod <- tmda(x[-test_index,], dims, y[-test_index],
          j, cov_str,
          maxit = iters, tol = 1e-10,
        )
        iteration <- mod$iteration
        pred <- predict.EM(x[test_index,], mod)$pred
        err <- loss(y[test_index], pred)
        inits[[countj]] <- list()
        inits[[countj]][["mu"]] <- mod$mu
        inits[[countj]][["pi"]] <- mod$pi
        inits[[countj]][["sig"]] <- mod$sig
        inits[[countj]][["gamma"]] <- mod$gamma
      }
      cv_error[countj, i] <- err
      iterations[countj, i] <- iteration
    }
  }
  cv_result <- list()
  cv_result[["cv_error"]] <- apply(cv_error, 1, mean)
  cv_result[["best subnum"]] <- subnum[which.min(cv_result[["cv_error"]])]
  cv_result[["iteration"]] <- iterations[,1]
  cv_result[["inits"]] <- inits
  print(paste("best number of mixtures: ", cv_result[["best subnum"]]))
  cv_result
}

getAIC = function(model, critera = "BIC"){
  n = sum(sapply(model$gamma, nrow))
  if(critera == "AIC"){
    k = 2
  }else if(critera == "BIC"){
    k = log(n)
  }else{
    stop("incorrect criterion")
  }
  subnum = sapply(model$pi, length)
  dims = dim(model$mu[[1]][[1]])
  if(model$cov_str == "diag"){
    p = sum(subnum) + prod(dims) * sum(subnum) + sum(dims)
  }else{
    p = sum(subnum) + prod(dims) * sum(subnum) + sum(dims * (dims + 1)/2)
  }
  return(-2 * model[["loglike"]] + k * p/n)
}