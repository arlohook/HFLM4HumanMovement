library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape)
library(fields)

sim.results <- readRDS("./Simulation Results/Cleaned/Pen Simulation Results 50.rds")
sim.data <- readRDS("Simulations Data 50.rds")

true.beta <- sim.data$Model$beta.est.list[[1]]
true.alpha <- sim.data$Model$intercept.fd

 # set up parameters for evaluating the surface
nfine <- 101
tfine <- seq(0, 1, length = nfine)
obspts <- cbind(rep(tfine,nfine),rep(tfine,1,each=nfine))

## get betas

bmat.list <- lapply(sim.results$Results, function(x){
  
  bm <- eval.FEM(FEM = x$beta.est.list[[1]], locations = obspts)
  
})



beta.mean <- bmat.list[[1]]

for (i in 2:1000) {
  
  beta.mean <- beta.mean + bmat.list[[i]]
  
}

beta.mean = matrix(beta.mean/1000, nfine, nfine)

true.b.mat <- matrix(eval.FEM(FEM = true.beta, obspts), nfine, nfine)
beta.res = true.b.mat - beta.mean
image.plot(beta.mean)

# get mean and sd SSIE
image.plot(beta.mean)
image.plot(beta.res)

# Function to perform Simpson's rule integration on the surface
SIE <- function(surface_matrix, h_x = 0.01, h_y = 0.01) {
  n <- nrow(surface_matrix)
  m <- ncol(surface_matrix)
  
  sq_matrix = surface_matrix^2
  sq_matrix[is.na(sq_matrix)] <- 0
  
  if (n %% 2 == 0 || m %% 2 == 0) {
    stop("Simpson's rule requires an odd number of intervals in both directions.")
  }
  
  # Simpson's rule coefficients
  simpson_coeff <- function(n) {
    coeff <- rep(1, n)
    coeff[seq(2, n-1, by=2)] <- 4
    coeff[seq(3, n-2, by=2)] <- 2
    return(coeff)
  }
  
  row_coeff <- simpson_coeff(n)
  col_coeff <- simpson_coeff(m)
  
  # Apply Simpson's rule
  integral <- sum(row_coeff %*% sq_matrix %*% col_coeff) * (h_x * h_y / 9)
  
  return(integral)
}


bSIE.list <- lapply(bmat.list, function(x){
  
  res <- true.b.mat - matrix(x, nfine, nfine)
  
  sie <- SIE(res)
  
  return(sie)
  
})

bSIE.vec <- unlist(bSIE.list)

bMSIE <- mean(bSIE.vec)
bSDSIE = sd(bSIE.vec)

ggplot()+
  geom_density(aes(x = bSIE.vec))

# get alpha

alpha.list <- lapply(sim.results$Results, function(x){
  x$intercept.fd
})

a.co <- c()

for (i in 1:1000) {
  
  a.co = cbind(a.co, alpha.list[[i]]$coefs[,1])
  
}

alpha.mean.fd <- mean.fd(fd(a.co, alpha.list[[1]]$basis))


plot(alpha.mean.fd)
plot(true.alpha)

aSIE.list <- lapply(alpha.list, function(x){
  
  s.res <- ((true.alpha-x)^2)
  cat(".")
  sie <- inprod(s.res)
  
  return(sie)
})




aSIE.vec <- unlist(aSIE.list)
aMSIE <- mean(aSIE.vec)
aSDSIE <- sd(aSIE.vec)

result.list = list("alpha" = alpha.mean.fd, "beta" = beta.mean,
                   "BMSE" = bMSIE, "BSDSE" = bSDSIE,
                   "AMSE" = aMSIE, "ASDSE" = aSDSIE)
saveRDS(result.list,file = "./Simulation Results/Cleaned/Pen 50 Processed.rds")

