## simulation 
library(refund)
library(fda)
library(tidyverse)
library(reshape2)
library(parallel)
library(foreach)
library(doRNG)
library(doParallel)
library(progressr)
library(progress)


setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")
# read in data set

Data = readRDS("./Simulation Study/Data/Simulations Data 500.rds")

n = ncol(Data$Data[[1]]$x.fd$coefs)

nfine = 101
tfine = seq(0, 1, length = nfine)


# function to extract coefficients

get.info <- function(model){
  
  coefs <- coef(model)
  alpha <- coefs$smterms$`Intercept(tfine)`$coef
  alpha$value = alpha$value + coefs$pterms[1]
  beta <- coefs$smterms$`ff(X)`$coef
  residuals <- model$residuals
  
  return(list("alpha" = alpha, "beta" = beta, "residuals" = residuals))
  
}

# function to fit the model

fit.sim = function(i){
  
  
  x.fd = eval.fd(tfine, Data[[3]][[i]]$x.fd)
  y.fd = eval.fd(tfine, Data[[3]][[i]]$yplusE)
  
  pffr.data = data.frame("ID" = 1:n)
  pffr.data$Y = t(y.fd)
  pffr.data$X = t(x.fd)
  
  
  A = Sys.time()
  
  fit = pffr(Y ~ ff(X, limits = "s<=t", xind = tfine,
                    splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                      k = c(10, 20))), 
             yind = tfine, data = pffr.data, bs.int = list(bs = 'ps', k = 20, m = c(3,2)))
  
  
  coef = get.info(fit)
  coef = append(coef, Sys.time()-A)
  names(coef)[length(coef)] = "Run Time"
  
  
  
  return(coef)
  
  
}



# set up cluster
n_cores <- 15
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(
  cl,
  c("Data", "n", "nfine", "tfine", "get.info", "fit.sim"),
  envir = environment()
)

# run simulation




  
  RESULTS  = foreach(L = 1:100,
                          .packages = c("refund", "mgcv", "fda"),
                          .combine = function(x, y) {

                              append(x, list(y))},
                     
                          .init = list()) %dopar% {
                            
                            
                            fit.sim(i = L)}
  
  

# stop the cluster

stopCluster(cl)

# save results


details = "pffr fit, k = c(10, 20), m = c(2, 1)"
sim_results <- list("Details" = details, "Results" = RESULTS)
saveRDS(sim_results, "./Simulation Study/Results/pffr Simulation Results 500.rds")

