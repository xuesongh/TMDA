subroutine aperm_f(tnsrOut, tnsr, dims, nlen, perm, ndim)
    !new implementation, faster than the old one in most of cases
    !Goal: permute dimension of a tensor, mainly used for mode-m matricization (order = [2:k,1,(k+1):ndim], this is not runable).
    !Input variable
        !tnsr: vectorized tensor, 
        !dims: permuted dimension of tensor, not original!!!, 
        !nlen: prod(dim(tensor)),
        !perm: permuted order, 
        !ndim: length(dim(tensor)).
    !Output variable: tnsr: vectorized tensor.
        Integer, intent(inout):: ndim
        Integer:: nlen
        Integer:: perm(ndim)
        Integer, intent(inout):: dims(ndim)
        Real*8:: tnsr(nlen)
        Real*8:: tnsrOut(nlen)
        Integer:: j
        Integer:: dimsr(ndim)
        Integer:: iip(ndim)
        Integer:: i
        Integer:: stride(ndim)
        Integer:: lj,li,itmp
        !--------- Start
        dimsr = dims(perm)
        iip = 1
        do i = 2,ndim
            iip(i) = iip(i-1)*dims(i-1)
        end do
        stride = iip(perm)
        iip = 0
        lj = 1
        do li = 1,nlen
            tnsrOut(li)=tnsr(lj)
            do itmp = 1,ndim
                if(iip(itmp) == dimsr(itmp)-1) then ! no -1
                    iip(itmp) = 0
                else
                    iip(itmp) = iip(itmp) + 1
                    exit
                end if
            end do
            lj = 1 + sum(iip * stride)
        end do
end subroutine aperm_f



subroutine aperm_c(tnsr, dims, nlen, perm, ndim)
    !new implementation, faster than the old one in most of cases
    !Goal: permute dimension of a tensor, mainly used for mode-m matricization (order = [2:k,1,(k+1):ndim], this is not runable).
    !Input variable
        !tnsr: vectorized tensor, 
        !dims: permuted dimension of tensor, not original!!!, 
        !nlen: prod(dim(tensor)),
        !perm: permuted order, 
        !ndim: length(dim(tensor)).
    !Output variable: tnsr: vectorized tensor.
        Integer, intent(inout):: ndim
        Integer:: nlen
        Integer:: perm(ndim)
        Integer, intent(inout):: dims(ndim)
        Real*8:: tnsr(nlen)
        Real*8:: tnsrOut(nlen)
        Integer:: j
        Integer:: dimsr(ndim)
        Integer:: iip(ndim)
        Integer:: i
        Integer:: stride(ndim)
        Integer:: lj,li,itmp
        !--------- Start
        dimsr = dims(perm)
        iip = 1
        do i = 2,ndim
            iip(i) = iip(i-1)*dims(i-1)
        end do
        stride = iip(perm)
        iip = 0
        lj = 1
        do li = 1,nlen
            tnsrOut(li)=tnsr(lj)
            do itmp = 1,ndim
                if(iip(itmp) == dimsr(itmp)-1) then ! no -1
                    iip(itmp) = 0
                else
                    iip(itmp) = iip(itmp) + 1
                    exit
                end if
            end do
            lj = 1 + sum(iip * stride)
        end do
        tnsr = tnsrOut
end subroutine aperm_c
        
subroutine inverse_f(tnsrOut, mat, newdims, nlen, ndim, k)
!Goal: reverse the mode-k matrix to original tensor
!Input variable
    !tnsrOut: vectorized tensor, 
    !mat: vectorized matrix, 
    !newdims: dimension of tensorOut (k-th dim is the first), 
    !nlen: prod(dim(tensor)),
    !ndim: length(dim(tnsr)),
    !k: mat is mode-k mat of old tensor
!Output variable: tnsrOut: vectorized tensor.
    implicit integer (d, n)
    implicit double precision (t, m)
    Integer, intent(inout):: ndim
    Integer:: nlen, k
    Integer, intent(inout):: newdims(ndim)
    Real*8:: mat(nlen)
    Real*8:: tnsrOut(nlen)
    ! temp var
    Integer:: perm(ndim)
    Integer:: olddims(ndim)
    ! start
    if(k == 1) then
        tnsrOut = mat
        return
    else if(k == ndims) then
        olddims = [newdims(k), newdims(1:(k-1))]
        perm = [(i, i = 2,k), 1]
    else
        olddims = [newdims(k), newdims(1:(k-1)), newdims((k+1):ndim)]
        perm = [(i, i = 2,k), 1, (i, i = k+1,ndim)]
    end if
    call aperm_f(tnsrOut, mat, olddims, nlen, perm, ndim)
    return
    end subroutine inverse_f
    
subroutine mode_k_mat_f(matOut, tnsr, olddims, nlen, ndim, k)
!Goal: the mode-k matrix of a tensor
!Input variable
    !tnsrOut: vectorized tensor, 
    !mat: vectorized matrix, 
    !olddims: dimension of tensorized mat (k-th dim is the first), 
    !nlen: prod(dim(tensor)),
    !ndim: length(dim(tnsr)),
    !k: mat is mode-k mat of old tensor
!Output variable: tnsrOut: vectorized tensor.
    implicit integer (d, n)
    implicit double precision (t, m)
    Integer, intent(inout):: ndim
    Integer:: nlen, k
    Integer, intent(inout):: olddims(ndim)
    Real*8:: tnsr(nlen)
    Real*8:: matOut(nlen)
    Integer:: i
    Integer:: neworder(ndim)
    ! start
    if(k == 1) then
        matOut = tnsr
        return
    else if(k == ndim) then
        neworder = [k, (i, i = 1,k-1)]
    else
        neworder = [k, (i, i = 1,k-1), (i, i = (k+1),ndim)]
    end if
    call aperm_f(matOut, tnsr, olddims, nlen, neworder, ndim)
    return
    end subroutine mode_k_mat_f
    
subroutine mode_k_prod_f(tnew, newdims, newlen, tnsr, tnsrdims, nlen, ndim, mat, matdims, k)
!Goal: mode k product between a tensor and a matrix.
!Input variable
    !tnew: vectorized tensor, used to store output,
    !newdims: dim(tnew),
    !newlen: prod(dim(tnew)),
    !tnsr: vectorized tensor, 
    !tnsrdims: dim(tnsr), 
    !nlen: prod(dim(tnsr)),
    !ndim: length(dim(tnsr)),
    !mat: matrix, 
    !matdims: dimension of matrix, 
    !k: mode k.
!Output variable: tnsr: vectorized tensor.
    implicit integer (d, n, k)
    implicit double precision (t, m)
    Integer:: k
    Integer:: ndim
    Integer:: nlen
    Integer:: newlen
    Integer:: tnsrdims(ndim)
    Integer:: matdims(2)
    Integer:: newdims(ndim)
    Real*8:: tnew(newlen)
    Real*8:: mat(matdims(1), matdims(2))
    Real*8:: tnsr(nlen)
    Integer:: i
    Real*8:: temp_tnsr(nlen)
    Integer:: neworder(ndim)
    ! start
    if(k == 1) then
        temp_tnsr = tnsr
    else if(k == ndim) then
        neworder = [k, (i, i = 1,k-1)]
        call aperm_f(temp_tnsr, tnsr, tnsrdims, nlen, neworder, ndim)
    else
        neworder = [k, (i, i = 1,k-1), (i, i = (k+1),ndim)]
        call aperm_f(temp_tnsr, tnsr, tnsrdims, nlen, neworder, ndim)
    end if
    call inverse_f(tnew, matmul(mat, reshape(temp_tnsr, [matdims(2), nlen / matdims(2)])),newdims, newlen, ndim, k)
    return
    end subroutine mode_k_prod_f

subroutine tucker_prod_f(tnew, newdims, tnsr, tnsrdims, ndim, longvec, veclen, skip)
    implicit none
    Integer:: ndim

    Integer:: veclen
    Integer:: tnsrdims(ndim)
    Integer:: newdims(ndim)
    Real*8:: tnew(product(newdims))
    Real*8, intent(inout):: longvec(veclen)
    Real*8, intent(inout):: tnsr(product(tnsrdims))
    Integer, optional:: skip(ndim)
    Integer:: lenIn
    Integer:: lenOut
    Integer:: matshape(2)
    Integer:: currentlen
    Integer:: dimIn(ndim)
    Integer:: dimOut(ndim)
    Integer:: pastlen
    Integer:: iter
    Integer:: nlen
    Integer:: newlen
    Integer:: io,i
    Integer:: maxlen
    Real*8:: tnsrIn(product((/(maxval((/tnsrdims(i), newdims(i)/)), i = 1,ndim)/)))
    Real*8:: tnsrOut(product((/(maxval((/tnsrdims(i), newdims(i)/)), i = 1,ndim)/)))
    Real*8:: currentmat(maxval(newdims)*maxval(tnsrdims))
    maxlen = size(tnsrIn)
    newlen = product(newdims)
    nlen = product(tnsrdims)
    pastlen = 0
    iter = 0
    lenIn = nlen
    dimIn = tnsrdims
    tnsrIn(:nlen) = tnsr 
    do while (iter < ndim)
        iter = iter + 1
        if(present(skip) .and. skip(iter) == 1) then
            if(mod(iter, 2) == 1) then
                dimOut = dimIn
                lenOut = lenIn / tnsrdims(iter) * newdims(iter)
                dimOut(iter) = newdims(iter)
                tnsrOut = tnsrIn
                io = 0   
            else
                dimIn = dimOut
                lenIn = lenOut / tnsrdims(iter) * newdims(iter)
                dimIn(iter) = newdims(iter)
                tnsrIn = tnsrOut
                io = 1
            end if
        pastlen = pastlen + product([newdims(iter), tnsrdims(iter)])
        cycle
        end if
        matshape = [newdims(iter), tnsrdims(iter)]
        currentlen = product(matshape)
        currentmat(:currentlen) = longvec((pastlen+1):(pastlen+currentlen))
        pastlen = pastlen + currentlen
        if(mod(iter, 2) == 1) then
            dimOut = dimIn
            lenOut = lenIn / matshape(2) * matshape(1)
            dimOut(iter) = matshape(1)
            call mode_k_prod_f(tnsrOut, dimOut, lenOut, tnsrIn, dimIn, lenIn, ndim, currentmat, matshape, iter)
            io = 0
        else
            dimIn = dimOut
            lenIn = lenOut / matshape(2) * matshape(1)
            dimIn(iter) = matshape(1)
            call mode_k_prod_f(tnsrIn, dimIn, lenIn, tnsrOut, dimOut, lenOut, ndim, currentmat, matshape, iter)
            io = 1
        end if
    end do
    if(io == 1) then
        tnew = tnsrIn
    else
        tnew = tnsrOut
    end if
    return
    end subroutine tucker_prod_f

subroutine mode_k_diag_prod_f(tnsr, tnsrdims, nlen, ndim, diag, k)
    !Goal: mode k product between a tensor and a **diagonal** matrix.
    !Input variable
        !tnsr: vectorized tensor, 
        !tnsrdims: dim(tnsr), 
        !nlen: prod(dim(tnsr)),
        !ndim: length(dim(tnsr)),
        !diag: vector, 
        !k: mode k.
    !Output variable: tnsr: vectorized tensor.
    implicit integer (d, n, k, p)
    implicit double precision (t, m)
    Integer:: k, j, d
    Integer:: ndim
    Integer:: nlen
    Integer:: tnsrdims(ndim)
    Real*8:: diag(tnsrdims(k))
    Real*8:: tnsr(nlen)
    Integer:: i
    Real*8:: temp_tnsr(nlen)
    Integer:: neworder(ndim)
    ! start
    if (ndim == 1) then
        do i = 1, nlen
            tnsr(i) = tnsr(i) * diag(i)
        end do
        return
    end if
    p_pri = product(tnsrdims(1:(k-1)))
    p_tok = p_pri*tnsrdims(k)
    p_post = product(tnsrdims((k+1):ndim))
    do d = 1, tnsrdims(k)
        tnsr([(((j-1)*p_tok+(d-1)*p_pri+i, i = 1,p_pri), j = 1,p_post)]) = &
        tnsr([(((j-1)*p_tok+(d-1)*p_pri+i, i = 1,p_pri), j = 1,p_post)]) * diag(d)
    end do
    return
    end subroutine mode_k_diag_prod_f

subroutine tucker_diag_prod_f(tnsr, tnsrdims, ndim, longvec, skip)
    implicit none
    Integer:: ndim
    Integer:: tnsrdims(ndim)
    Real*8, intent(inout):: longvec(sum(tnsrdims))
    Real*8:: tnsr(product(tnsrdims))
    Integer, optional:: skip(ndim)
    Real*8:: currentdiag(maxval(tnsrdims))
    Integer:: currentlen
    Integer:: pastlen
    Integer:: iter
    Integer:: nlen
    nlen = product(tnsrdims)
    pastlen = 0
    iter = 0
    do while (iter < ndim)
        iter = iter + 1
        if(present(skip) .and. skip(iter) == 1) then
        pastlen = pastlen + tnsrdims(iter)
        cycle
        end if
        currentlen = tnsrdims(iter)
        currentdiag(:currentlen) = longvec((pastlen+1):(pastlen+currentlen))
        pastlen = pastlen + currentlen
        call mode_k_diag_prod_f(tnsr, tnsrdims, nlen, ndim, currentdiag, iter)
    end do
    end subroutine tucker_diag_prod_f

subroutine EM_diag(longx, xdims, ndim, y, n, nclass, subnum, longmu, diagsig_inv_sq, &
                    pis, longgamma, glen, maxiter, loglike)
    implicit none
    interface
        subroutine tucker_diag_prod_f(tnsr, tnsrdims, ndim, longvec, skip)
            Integer:: ndim
            Integer:: tnsrdims(ndim)
            Real*8, intent(inout):: longvec(sum(tnsrdims))
            Real*8:: tnsr(product(tnsrdims))
            Integer, optional:: skip(ndim)
        end subroutine
    end interface
    Integer:: maxiter
    Integer:: ndim
    Integer:: n
    Integer:: nclass
    Integer:: subnum(nclass)
    Integer:: y(n)
    Integer:: xdims(ndim)
    Integer:: glen
    Real*8:: longx(product([xdims, n]))
    Real*8:: longmu(product([xdims, sum(subnum)]))
    Real*8:: pis(sum(subnum))
    Real*8:: longgamma(glen)
    Real*8:: loglike(maxiter)
    Real*8:: diagsig_inv_sq(sum(xdims))
    Integer:: iter, j, i, r, pastx, m
    Integer:: pastsig
    Real*8:: dist(maxval(subnum))
    Real*8:: tempIn(product(xdims))
    Real*8:: tempOut(product(xdims))
    Integer:: xlen
    Integer:: skip(ndim)
    Real*8:: U(maxval(xdims))
    Real*8:: temp_mat(product(xdims))
    Integer:: currentg, pastg, currentmu, pastmu, currentpi, pastpi
    Real*8:: logdets(ndim)
    Real*8:: mindist
    !start
    xlen = size(tempIn)
    logdets = 0.
    do iter = 1, maxiter
        loglike(iter) = 0.
        if(mod(iter, (maxiter+4)) == 0) print*, "EM_diag ", "iter: ", iter
        pastg = 0
        pastx = 0
        !E-step: update gamma
        do i = 1, n
            currentg = subnum(y(i))
            pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
            pastpi = sum(subnum(:y(i))) - subnum(y(i))
            do r = 1, currentg
                tempOut = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                pastmu = pastmu + xlen
                call tucker_diag_prod_f(tempOut, xdims, ndim, diagsig_inv_sq, skip)
                dist(r) = dot_product(tempOut, tempOut)
            end do
            pastx = pastx + xlen
            mindist = minval(dist(:currentg))
            dist(:currentg) = dist(:currentg) - minval(dist(:currentg))
            dist(:currentg) = exp(-dist(:currentg) / 2) * pis((pastpi+1) : (pastpi+subnum(y(i)))) + 1e-15
            loglike(iter) = loglike(iter) + log(sum(dist(:currentg))+ 1e-15)-mindist
            longgamma((pastg + 1):(pastg+currentg)) = dist(:currentg) / sum(dist(:currentg))
            pastg = pastg + currentg
        end do
        if(any(isnan(longgamma))) then
            print*, "gamma nan"
            exit
        end if
        !M-step
        ! update pi and mu
        pis = 0.
        longmu = 0.
        pastg = 0
        pastx = 0
        do i = 1, n
            pastpi = sum(subnum(:y(i))) - subnum(y(i))
            currentg = subnum(y(i))
            pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
            do r = 1, currentg
                longmu((pastmu + 1) : (pastmu + xlen)) = longmu((pastmu + 1) : (pastmu + xlen)) +&
                    longx((pastx+1):(pastx+xlen)) * longgamma(pastg + r)
                pastmu = pastmu + xlen
            end do
            pastx = pastx + xlen
            pis((pastpi+1) : (pastpi+subnum(y(i)))) = pis((pastpi+1) : (pastpi+subnum(y(i)))) + &
                longgamma((pastg + 1):(pastg+currentg))
            pastg = pastg + currentg
        end do
        pastmu = 0
        do j = 1, nclass
            pastpi = sum(subnum(:j)) - subnum(j)
            do r = 1, subnum(j)
                longmu((pastmu + 1) : (pastmu + xlen)) = longmu((pastmu + 1) : (pastmu + xlen)) /&
                        pis(pastpi+r)
                pastmu = pastmu + xlen
            end do
            pis((pastpi+1):(pastpi+subnum(j))) = pis((pastpi+1):(pastpi+subnum(j))) /&
                sum(pis((pastpi+1):(pastpi+subnum(j))))
        end do
        if(any(isnan(pis))) then
            print*, "pis nan1"
            exit
        end if
        if(any(isnan(longmu))) then
            print*, "mu nan1"
            exit
        end if
        !update sigma
        pastsig = 0
        do m = 1, ndim
            skip = 0
            skip(m) = 1
            pastg = 0
            pastx = 0
            U = 0.
            do i = 1, n
                currentg = subnum(y(i))
                pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
                pastpi = sum(subnum(:y(i))) - subnum(y(i))
                do r = 1, currentg
                    tempOut = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                    pastmu = pastmu + xlen
                    call tucker_diag_prod_f(tempOut, xdims, ndim, diagsig_inv_sq, skip)
                    call mode_k_mat_f(temp_mat, tempOut, xdims, xlen, ndim, m)
                    U(:xdims(m)) = U(:xdims(m)) + &
                    norm2(reshape(temp_mat, [xdims(m), xlen / xdims(m)]), dim=2) ** 2 * longgamma(pastg+r)
                end do
                pastx = pastx + xlen
                pastg = pastg + currentg
            end do  
            if (m < ndim) then
                diagsig_inv_sq((pastsig + 1): (pastsig + xdims(m))) = U(:xdims(m)) / maxval(U(:xdims(m)))
            else
                diagsig_inv_sq((pastsig + 1): (pastsig + xdims(m))) = U(:xdims(m)) / (n * xlen / xdims(m))
            end if
            logdets(m) = log(product(diagsig_inv_sq((pastsig + 1): (pastsig + xdims(m))))+1e-15) * (-xlen / xdims(m) / 2)
            diagsig_inv_sq((pastsig + 1): (pastsig + xdims(m))) = 1 / sqrt(diagsig_inv_sq((pastsig+1):(pastsig+xdims(m)))+1e-15)
            pastsig = pastsig + xdims(m)
        end do
        if(any(isnan(diagsig_inv_sq))) then
            print*, "sig nan"
            exit
        end if
        loglike(iter) = loglike(iter)/n + sum(logdets)
        if ( (iter > 1) .and. (iter < maxiter) .and. (abs((loglike(iter) - loglike(iter - 1))&
            / loglike(iter - 1)) < loglike(maxiter))) then
            maxiter = iter
            print*, "EM_diag ", "early stop at iter: ", iter
            exit
        end if
    end do
    end subroutine EM_diag

subroutine pred_diag(longx, xdims, ndim, prob, n, nclass, subnum, longmu, diagsig_inv_sq, pis, phi)
    implicit none
    interface
        subroutine tucker_diag_prod_f(tnsr, tnsrdims, ndim, longvec, skip)
            Integer:: ndim
            Integer:: tnsrdims(ndim)
            Real*8, intent(inout):: longvec(sum(tnsrdims))
            Real*8:: tnsr(product(tnsrdims))
            Integer, optional:: skip(ndim)
        end subroutine
    end interface
    Integer:: ndim
    Integer:: n
    Integer:: nclass
    Integer:: subnum(nclass)
    Real*8:: prob(n, nclass)
    Integer:: xdims(ndim)
    Real*8:: longx(product([xdims, n]))
    Real*8:: longmu(product([xdims, sum(subnum)]))
    Real*8:: pis(sum(subnum))
    Real*8:: phi(nclass)
    Real*8:: diagsig_inv_sq(sum(xdims))
    Integer:: iter, j, i, r, pastx, currentmu, pastmu, currentpi, pastpi, m 
    Real*8:: dist(sum(subnum))
    Real*8:: tempOut(product(xdims))
    Integer:: xlen
    !start
    xlen = size(tempOut)
    pastx = 0
    do i = 1, n
        pastmu = 0
        pastpi = 0
        do j = 1,nclass 
            do r = 1, subnum(j)
                tempOut = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                pastmu = pastmu + xlen
                call tucker_diag_prod_f(tempOut, xdims, ndim, diagsig_inv_sq)
                dist(pastpi + r) = sum(tempOut * tempOut)
            end do
            pastpi = pastpi + subnum(j)
        end do
        dist = exp(- (dist - minval(dist)) / 2) * pis
        pastpi = 0
        do j = 1,nclass
            prob(i, j) = phi(j) * sum(dist((pastpi+1):(pastpi+subnum(j))))
            pastpi = pastpi + subnum(j)
        end do
        prob(i,:) = prob(i,:) / sum(prob(i,:))
        pastx = pastx + xlen
    end do
end subroutine pred_diag

subroutine EM(longx, xdims, ndim, y, n, nclass, subnum, longmu, longsig_inv_sq, pis, longgamma, glen, maxiter, loglike)
    implicit none
    interface
        subroutine tucker_prod_f(tnew, newdims, tnsr, tnsrdims, ndim, longvec, veclen, skip)
            Integer:: ndim
            Integer:: veclen
            Integer:: tnsrdims(ndim)
            Integer:: newdims(ndim)
            Real*8:: tnew(product(newdims))
            Real*8, intent(inout):: longvec(veclen)
            Real*8, intent(inout):: tnsr(product(tnsrdims))
            Integer, optional:: skip(ndim)
        end subroutine
    end interface
    Integer:: maxiter
    Integer:: ndim
    Integer:: n
    Integer:: nclass
    Integer:: subnum(nclass)
    Integer:: y(n)
    Integer:: xdims(ndim)
    Integer:: glen
    Real*8:: longx(product([xdims, n]))
    Real*8:: longmu(product([xdims, sum(subnum)]))
    Real*8:: pis(sum(subnum))
    Real*8:: longgamma(glen)
    Real*8:: loglike(maxiter)
    Real*8:: longsig_inv_sq(sum(xdims ** 2))
    Integer:: iter, j, currentg, pastg, i, r, pastx, currentmu, pastmu, currentpi, pastpi, m
    Integer:: pastsig
    Real*8:: dist(maxval(subnum))
    Real*8:: tempIn(product(xdims))
    Real*8:: tempOut(product(xdims))
    Integer:: xlen, lworks, info
    Integer:: skip(ndim)
    Real*8:: S_mat(maxval(xdims),maxval(xdims))
    Real*8:: S(maxval(xdims))
    Real*8:: U(maxval(xdims),maxval(xdims))
    Real*8, allocatable:: temp_mat(:,:)
    Real*8, allocatable:: works(:)
    Real*8:: logdets(ndim)
    Real*8:: mindist
    xlen = size(tempIn)
    logdets = 0.
    !start
    do iter = 1, maxiter
        if(mod(iter, (maxiter+4)) == 0) print*, "EM_general ", "iter: ", iter
        loglike(iter) = 0.
        pastg = 0
        pastx = 0
        !E-step: update gamma
        do i = 1, n
            currentg = subnum(y(i))
            pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
            pastpi = sum(subnum(:y(i))) - subnum(y(i))
            do r = 1, currentg
                tempIn = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                pastmu = pastmu + xlen
                call tucker_prod_f(tempOut, xdims, tempIn, xdims, ndim, longsig_inv_sq, size(longsig_inv_sq))
                dist(r) = dot_product(tempOut, tempOut)
            end do
            pastx = pastx + xlen
            mindist = minval(dist(:currentg))
            dist(:currentg) = dist(:currentg) - minval(dist(:currentg))
            dist(:currentg) = exp(-dist(:currentg) / 2) * pis((pastpi+1) : (pastpi+subnum(y(i))))
            loglike(iter) = loglike(iter) + log(sum(dist(:currentg))+ 1e-15)-mindist
            longgamma((pastg + 1):(pastg+currentg)) = dist(:currentg) / sum(dist(:currentg))
            pastg = pastg + currentg
        end do
        !M-step
        ! update pi and mu
        pis = 0.
        longmu = 0.
        pastg = 0
        pastx = 0
        pastsig = 0
        do i = 1, n
            pastpi = sum(subnum(:y(i))) - subnum(y(i))
            currentg = subnum(y(i))
            pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
            do r = 1, currentg
                longmu((pastmu + 1) : (pastmu + xlen)) = longmu((pastmu + 1) : (pastmu + xlen)) +&
                    longx((pastx+1):(pastx+xlen)) * longgamma(pastg + r)
                pastmu = pastmu + xlen
            end do
            pastx = pastx + xlen
            pis((pastpi+1) : (pastpi+subnum(y(i)))) = pis((pastpi+1) : (pastpi+subnum(y(i)))) + &
                longgamma((pastg + 1):(pastg+currentg))
            pastg = pastg + currentg
        end do
        pastmu = 0
        do j = 1, nclass
            pastpi = sum(subnum(:j)) - subnum(j)
            do r = 1, subnum(j)
                longmu((pastmu + 1) : (pastmu + xlen)) = longmu((pastmu + 1) : (pastmu + xlen)) /&
                        pis(pastpi+r)
                pastmu = pastmu + xlen
            end do
            pis((pastpi+1):(pastpi+subnum(j))) = pis((pastpi+1):(pastpi+subnum(j))) /&
                sum(pis((pastpi+1):(pastpi+subnum(j))))
        end do
        !update sigma
        do m = 1, ndim
            skip = 0
            skip(m) = 1
            pastg = 0
            pastx = 0
            lworks = 3 * xdims(m)
            allocate(temp_mat(xdims(m), xlen / xdims(m)))
            allocate(works(lworks))
            U = 0.
            do i = 1, n
                currentg = subnum(y(i))
                pastmu = xlen * ( sum(subnum(:y(i))) - subnum(y(i)) )
                pastpi = sum(subnum(:y(i))) - subnum(y(i))
                
                do r = 1, currentg
                    tempIn = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                    pastmu = pastmu + xlen
                    call tucker_prod_f(tempOut, xdims, tempIn, xdims, ndim, longsig_inv_sq, size(longsig_inv_sq), skip)
                    call mode_k_mat_f(temp_mat, tempOut, xdims, xlen, ndim, m) 
                    U(:xdims(m), :xdims(m)) = U(:xdims(m), :xdims(m)) + &
                                            matmul(temp_mat, &
                                            transpose(temp_mat)) * &
                                            longgamma(pastg+r)   
                end do
                pastx = pastx + xlen
                pastg = pastg + currentg
            end do
            if (m < ndim) then
                longsig_inv_sq((pastsig + 1) : (pastsig + xdims(m)**2)) = [U(:xdims(m), :xdims(m))] / U(1,1)
            else
                longsig_inv_sq((pastsig + 1) : (pastsig + xdims(m)**2)) = [U(:xdims(m), :xdims(m))] / (n * xlen / xdims(m))
            end if
            U(:xdims(m), :xdims(m)) = reshape(longsig_inv_sq((pastsig + 1) : (pastsig + xdims(m)**2)), [xdims(m), xdims(m)])
            call dsyev("V", "U", xdims(m), U(:xdims(m), :xdims(m)), xdims(m), S(:xdims(m)), works, lworks, info)
            logdets(m) = log(product(S(:xdims(m)))+1e-15) * (-xlen / (2*xdims(m)))
            S_mat = 0.
            do i = 1, xdims(m)
                S_mat(:xdims(m),i) = U(:xdims(m),i) * 1/sqrt(S(i) + 1e-10)
            end do 
            longsig_inv_sq((pastsig + 1) : (pastsig + xdims(m)**2)) = &
                                            [matmul(S_mat(:xdims(m), :xdims(m)), transpose(U(:xdims(m), :xdims(m))))] 
            pastsig = pastsig + xdims(m)**2
            deallocate(temp_mat)
            deallocate(works)
        end do
        loglike(iter) = loglike(iter)/n + sum(logdets)
        if ( (iter > 1) .and. (iter < maxiter) .and. (abs(loglike(iter) - loglike(iter - 1) &
            / loglike(iter - 1)) < loglike(maxiter))) then
            maxiter = iter
            print*, "EM_general ", "early stop at iter: ", maxiter
            exit
        end if   
    end do
    end subroutine EM

subroutine pred_f(longx, xdims, ndim, prob, n, nclass, subnum, longmu, longsig_inv_sq, pis, phi)
    implicit none
    interface
        subroutine tucker_prod_f(tnew, newdims, tnsr, tnsrdims, ndim, longvec, veclen, skip)
            Integer:: ndim
            Integer:: veclen
            Integer:: tnsrdims(ndim)
            Integer:: newdims(ndim)
            Real*8:: tnew(product(newdims))
            Real*8, intent(inout):: longvec(veclen)
            Real*8, intent(inout):: tnsr(product(tnsrdims))
            Integer, optional:: skip(ndim)
        end subroutine
    end interface
    Integer:: ndim
    Integer:: n
    Integer:: nclass
    Integer:: subnum(nclass)
    Real*8:: prob(n, nclass)
    Integer:: xdims(ndim)
    Real*8:: longx(product([xdims, n]))
    Real*8:: longmu(product([xdims, sum(subnum)]))
    Real*8:: pis(sum(subnum))
    Real*8:: phi(nclass)
    Real*8:: longsig_inv_sq(sum(xdims ** 2))
    Integer:: iter, j, i, r, pastx, currentmu, pastmu, currentpi, pastpi, m 
    Real*8:: dist(sum(subnum))
    Real*8:: tempIn(product(xdims))
    Real*8:: tempOut(product(xdims))
    Integer:: xlen
    !start
    xlen = size(tempIn)
    pastx = 0
    do i = 1, n
        pastmu = 0
        pastpi = 0
        do j = 1,nclass 
            do r = 1, subnum(j)
                tempIn = longx((pastx+1):(pastx+xlen)) - longmu((pastmu + 1) : (pastmu + xlen))
                pastmu = pastmu + xlen
                call tucker_prod_f(tempOut, xdims, tempIn, xdims, ndim, longsig_inv_sq, size(longsig_inv_sq))
                dist(pastpi + r) = dot_product(tempOut, tempOut)
            end do
            pastpi = pastpi + subnum(j)
        end do
        dist = exp(- (dist - minval(dist)) / 2) * pis
        pastpi = 0
        do j = 1,nclass
            prob(i, j) = phi(j) * sum(dist((pastpi+1):(pastpi+subnum(j))))
            pastpi = pastpi + subnum(j)
        end do
        prob(i,:) = prob(i,:) / sum(prob(i,:)) 
        pastx = pastx + xlen
    end do
    end subroutine pred_f