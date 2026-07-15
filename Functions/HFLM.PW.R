

HFLM.PW = function(y.fd, x.fd, tn = 101, Ak = 20, Bk = c(20,20)){
  
  
  yr = y.fd$basis$rangeval
  
  t = seq(yr[1],yr[2], length = tn)
  
  Y = t(eval.fd(t, y.fd))
  
  X1 = t(eval.fd(t[1], x.fd))
  
  A = numeric(tn)
  B = matrix(NA, tn, tn)
  
  zmod = lm(Y[,1] ~ X1)
  
  A[1] = coef(zmod)[1]
  B[1,1] = coef(zmod)[2]
  
  kops = round(seq(5, Bk[1], length = 4),0)
  
  for(i in 2:tn){
    
    Yt = Y[,i]
    
    s = seq(0,t[i], length = tn)
    
    Xs = t(eval.fd(s, x.fd))
    
    if(i/tn <= 0.1){
      m = 0
      k = 2
      }else if(i/tn > 0.1 & i/tn <= 0.25){
          k = kops[1]
          m = 2
          }else if(i/tn > 0.25 & i/tn <= 0.5){
            k = kops[2]
            m = 2
            }else if(i/tn > 0.5 & i/tn <= 0.75){
              k = kops[3]
              m = 2
            }else{
                k = kops[4]
                m = 2
              }
    
    mod = pfr(Yt ~ lf(Xs, argvals = s, k = k, m = m, bs = 'ps'))
   
    co = coef(mod)
    
    A[i] = mod$coefficients[1]
    
    B[1:i,i] = approx(x = co$Xs.argvals, y = co$value, xout = t[1:i])$y
    
  }
  
  
  A2sm = data.frame("A" = A, t = t)
  B2sm = data.frame("B" = melt(B)$value, "s" = rep(t, tn), "t" = rep(t, each = tn)) %>% na.omit()
  
  
  Asm.mod = bam(A ~ s(t, k = Ak, m = c(2,1)), data = A2sm)
  Bsm.mod = bam(B ~ te(s, t, k = Bk, bs = c("ps", "ps"), m = c(2,1)), data = B2sm)
  
  Bsm = data.frame("B" = predict(Bsm.mod, newdata = B2sm), "s" = B2sm$s, "t" = B2sm$t)
  Asm = data.frame("A" = predict(Asm.mod, newdata = A2sm), "t" = A2sm$t)
  
  out = list("alpha" = Asm, "beta" = Bsm, "beta.raw" = B2sm, "alpha.raw" = A2sm)
  
  return(out)
  
  
}