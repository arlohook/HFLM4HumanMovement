
FEM.diff <- function(grid, nderiv, axis){
  
  if(axis == "s"){
    
  
    if(nderiv == 1){
      
        nfine = nrow(grid)
      
        A = (grid[2:nfine,] - grid[1:(nfine-1),])/(1/nfine)
      
        X = A[,-1]
        b = diag(X)
      
        d = A[nrow(A), ncol(A)]+(A[nrow(A), ncol(A)] - A[(nrow(A)-1), ncol(A)])
      
        e = c(rep(NA, nrow(A)))
        f = c(b,d)
        B = rbind(A, e)
      
      
        diag(B) <- f
        
        B[is.na(B)] = 0
        
      return(B)
      }
    
    if(nderiv == 2){
      
      # 1st
      nfine = nrow(grid)
      
      A = (grid[2:nfine,] - grid[1:(nfine-1),])/(1/nfine)
      
      X = A[,-1]
      b = diag(X)
      
      d = A[nrow(A), ncol(A)]+(A[nrow(A), ncol(A)] - A[(nrow(A)-1), ncol(A)])
      
      e = c(rep(NA, nrow(A)))
      f = c(b,d)
      B = rbind(A, e)
      
      
      diag(B) <- f
      
      # 2nd
      nfine = nrow(B)
      
      A2 = (B[2:nfine,] - B[1:(nfine-1),])/(1/nfine)
      
      X2 = A2[,-1]
      b2 = diag(X2)
      
      d2 = A[nrow(A2), ncol(A2)]+(A2[nrow(A2), ncol(A2)] - A2[(nrow(A2)-1), ncol(A2)])
      
      e2 = c(rep(NA, nrow(A2)))
      f2 = c(b2,d2)
      B2 = rbind(A2, e2)
      
      
      diag(B2) <- f2
      
      B2[is.na(B2)] = 0
      
      return(B2)
      
    }
    
    }
  
  
  
  
  if(axis == "t"){
    
    
    if(nderiv == 1){
      
      nfine = nrow(grid)
      
      A = (grid[,2:nfine] - grid[,1:(nfine-1)])/(1/nfine)
      
      
      X = A[-nrow(A),]
      
      b = diag(X)
      
      d = X[1,1]
      
      e = c(rep(NA, nrow(A)))
      f = c(d,b)
      B = cbind(e, A)
      
      
      diag(B) <- f
      
      B[is.na(B)] = 0
      
      return(B)
    }
    
    if(nderiv == 2){
      
      # 1st
      
      nfine = nrow(grid)
      
      A = (grid[,2:nfine] - grid[,1:(nfine-1)])/(1/nfine)
      
      
      X = A[-nrow(A),]
      
      b = diag(X)
      
      d = X[1,1]
      
      e = c(rep(NA, nrow(A)))
      f = c(d,b)
      B = cbind(e, A)
      
      
      diag(B) <- f
      
      # 2nd
      
      nfine = nrow(B)
      
      A2 = (B[,2:nfine] - B[,1:(nfine-1)])/(1/nfine)
      
      
      X2 = A2[-nrow(A2),]
      
      b2 = diag(X)
      
      d2 = X2[1,1]
      
      e2 = c(rep(NA, nrow(A2)))
      f2 = c(d2,b2)
      B2 = cbind(e2, A2)
      
      
      diag(B2) <- f2
      
      B2[is.na(B2)] = 0
      
      return(B2)
      
    }
  }
}







