HFLM.smidge <- function(y.fd, xfdlist, betalist, nfine = 101, lambda, fit.int = T, idtol = 1e6, fem.lfd = 2) {
  
  
  
  if(is.null(nfine)){
    nfine = betalist[[2]]$FEMbasis$nbasis*length(betalist)
  }
  if (!is.fd(y.fd)) {
    stop("y.fd is not a functional data object.")
  }
  
  if(fit.int){
    int.basis = betalist[[1]]
    
    betalist = betalist[2:length(betalist)]
  }
  
  if(length(xfdlist) != length(betalist)){
    stop("number of basis objects does not match number of predictors")
  }
  
  if((nfine%%2) == 0){
    stop("nfine must be an odd number")
  }
  
  
  if(nfine < (1.5*y.fd$basis$nbasis)){
    warning("nfine should be greater than 1.5 times nbasis for yfd")
  }
  
  yran = y.fd$basis$rangeval
  
  tfine = seq(yran[1], yran[2], length = nfine)
  
  # build up Z matrix for linear regression
  
  
  
  
  xpars = length(xfdlist)
  
  ZZ = c()
  
  for(p in 1:xpars) {
    
    
    x.fd = xfdlist[[p]]
    
    if (!is.fd(x.fd)) {
      stop("x.fd is not a functional data object.")
    }
    
    xran = x.fd$basis$rangeval
    
    if(sum(xran-yran) != 0){
      stop("x.fd and y.fd ranges do not match")
    }
    
    # evaluate X over tfine
    
    Xvals = eval.fd(tfine, x.fd)
    
    basis = betalist[[p]]
    
    
    Z = c()
    
    
    
    for(i in 1:nfine){
      
      
      # Evaluation times at t = tfine[i]
      int_pts =  cbind(tfine, tfine[i])
      
      FEMvals = eval.FEM(basis, int_pts) %>% replace_na(0)
      
      
      # We want to form Z_ij  =  sum x_i[k] phi_j[k] ~ \int x_i(t)phi_j(t)dt
      # integral is numerically approximated using Simpson's Rule
      
      
      s = 1/3 * ((yran[2]-yran[1])/(nfine-1)) * c(1,rep(c(4,2),(nfine-3)/2),4,1)
      z = t(Xvals)%*%diag(s)%*%FEMvals
      
      
      Z = rbind(Z,z)
      
    }
    
    ZZ = cbind(ZZ,Z)
  }
  
  
  # make Y vals into a long vector
  
  Yvals = eval.fd(tfine, y.fd)
  
  YY = c(t(Yvals))
  
  n = ncol(Yvals)
  
  # add in intercept
  if(fit.int){
    
    int.eval = eval.basis(tfine, int.basis)
    
    int.cols = apply(int.eval, 2, function(J){sapply(J, function(N){rep(N, n)})})
    
    ZZ = cbind(int.cols, ZZ)
  }
  
  nb = sum(unlist(lapply(betalist, function(x){x$FEMbasis$nbasis})))
  
  if(sum(lambda) == 0 ){                                       
    
    bigPmat = matrix(0,ncol(ZZ), ncol(ZZ))
    
    XTX = t(ZZ)%*%ZZ
    
    eigs = eigen(XTX)$values
    cond.no = (max(eigs)/min(eigs))
    eig.check = cond.no < idtol
    
    if(eig.check){
      
      coefs = solve(XTX,  t(ZZ)%*%YY)
      
    }else{stop("Error: Model is unidentifiable - try reducing nbasis or increasing lambda.")}
    
    
  }else{
    
    # penalised version
    
    obspts = cbind(rep(tfine, nfine),rep(tfine, 1, each=nfine)) 
    b.eval = eval.FEM(FEM = basis, obspts)
    Pmat = c()
    
    for (p in 1:xpars) {
      
      FEMvals.ds = apply(b.eval, 2, function(x){
        
        node = matrix(x, nfine, nfine, byrow = F) %>% replace_na(0)
        
        d.node = FEM.diff(node, nderiv = fem.lfd, axis = "s")
        
        return(c(d.node))
      })
      
      FEMvals.dt = apply(b.eval, 2, function(x){
        
        node = matrix(x, nfine, nfine, byrow = F) %>% replace_na(0)
        
        d.node = FEM.diff(node, nderiv = fem.lfd, axis = "t")
        
        return(c(d.node))
        
      })
      
      P = (t(FEMvals.ds)%*%FEMvals.ds + t(FEMvals.dt)%*%FEMvals.dt)/(nfine-1)
      
      
      Zeros = matrix(0, nb/xpars, nb/xpars)
      
      
      
      H = matrix(c(rep(Zeros, p-1)), nb/xpars, (p-1)*(nb/xpars))
      
      G = matrix(c(rep(Zeros, xpars-p)), nb/xpars, (xpars-p)*(nb/xpars))
      
      P1 = cbind(H, P, G)
      
      Pmat = rbind(Pmat, P1)
      
    }
    
    if(fit.int){
      
      d.int.basis = eval.basis(tfine, int.basis, Lfdobj = 2)
      int.P = t(d.int.basis)%*% d.int.basis
      inb = int.basis$nbasis
      
      bigPmat = rbind(  matrix(0, inb, inb + nb),
                        cbind(matrix(0, nb, inb), Pmat) )
      
      bigPmat[1:inb, 1:inb] = int.P
      
    }else{bigPmat = Pmat}
    
    
    lvec = c(rep(lambda[1],inb), rep(lambda[2], nb))
    
    # correct ridge penalty for quadratic elements
    
    obspts = cbind(rep(tfine, nfine), rep(tfine, each = nfine))
    
    s = 1/3 * ((yran[2]-yran[1])/(nfine-1)) * c(1,rep(c(4,2),(nfine-3)/2),4,1)
    
    rp = apply(eval.FEM(FEM = FEMbasis, locations = obspts), 2, function(b){
      
      K = matrix(b, nfine, nfine)
      K[is.na(K)] = 0
      K = K^2
      Ks = rowSums(diag(s)%*%K)
      Kt = sum(s*Ks)
      
      return(Kt)
      
    })
    
    ridge.vec = c(rep(0,inb), rp*lambda[3])
    XTXP =  t(ZZ)%*%ZZ + diag(lvec)*bigPmat + diag(ridge.vec)
    
    eigs = eigen(XTXP)$values
    cond.no = max(eigs)/min(eigs)
    eig.check = (max(eigs)/min(eigs)) < idtol
    
    if(eig.check){
      
      coefs = solve(XTXP,  t(ZZ)%*%YY)
    }else{stop("Error: Model is unidentifiable - try reducing nbasis or increasing lambda.")}  
    
  }
  
  
  # extract beta coefs and make FEM fd
  if(fit.int){
    inb = int.basis$nbasis
    k = 1:inb
    
    beta = matrix(c(coefs[-k]), (nb/xpars), xpars)
  }else{
    beta = matrix(c(coefs), (nb/xpars), xpars)
  }
  Beta.FD.List <- vector("list", length = xpars)
  
  for (j in 1:xpars) {
    
    basis = betalist[[j]]
    
    b = beta[,j]
    
    Beta.FD.List[[j]] <- FEM(coeff = b, FEMbasis = basis$FEMbasis)
    
    
    
  }
  
  
  # extract intercept and make fd
  if(fit.int){
    int_co = coefs[1:inb]
    
    Int.fd = fd(int_co, int.basis)
    
  }else{Int.fd = NULL}
  
  
  # get yhat values
  
  yhat <- matrix(c(ZZ %*% coefs), nfine, n, byrow = T)
  
  Yhat_exp <- smooth.basis(argvals = tfine, y = yhat, fdParobj = y.fd$basis)
  
  Yhat.fd <- Yhat_exp$fd
  
  
  Model <- list(y.fd, xfdlist, betalist, Beta.FD.List, Int.fd, Yhat.fd, lambda, nfine, cond.no)
  
  names(Model) <- c("y.fd", "xfdlist", "betalist", "beta.est.list", "intercept.fd", "yhat.fd", "lambda", "nfine", "cond.no")
  
  
  return(Model)
  
}
