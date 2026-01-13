###############################################################################
# predict.FEM.HFLM
###############################################################################


predict.FEM.HFLM <- function(model, xfdlist, nfine) {
 

if(length(xfdlist) != length(model$xfdlist)){
  stop("xfdlist is the incorrect length")
}

yran = model$intercept.fd$basis$rangeval

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
  
  basis = model$betalist[[p]]
  
  
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

# add in intercept

n = ncol(xfdlist[[1]]$coefs)

ZZ = cbind(diag(length(tfine))%x%matrix(1,n,1), ZZ)


## get vector of coefficent values
  b = c(eval.fd(tfine, model$intercept.fd))
  
  for (i in 1:xpars) {
    
   k = (model$beta.est.list[[i]]$coeff)
   
   b = append(b, k)
    
  }

## Make predictions
  

  
  yhat <- matrix(c(ZZ %*% b), nfine, n, byrow = T)
  
  
  if(nrow(yhat) > model$y.fd$basis$nbasis){
    
    yhat.fd <- smooth.basis(argvals = tfine, y = yhat, fdParobj = model$y.fd$basis)$fd
  }else{yhat.fd = fd(yhat, y.fd$basis)}
  
  
  
  predictions <- list(yhat.fd, xfdlist, model)
  
  names(predictions) <- c("yhat.fd", "xfdlist", "model")


return(predictions)
  
}