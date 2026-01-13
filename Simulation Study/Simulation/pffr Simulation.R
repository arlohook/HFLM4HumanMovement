## simulation 
library(refund)
library(fda)
library(tidyverse)
library(reshape2)



# read in data set

Data = readRDS("Simulations Data 100.rds")

nfine = 101
tfine = seq(0, 1, length = nfine)

Fit.List <- vector('list', 1000)

get.info <- function(model){
  
  coefs <- coef(model)
  alpha <- coefs$smterms$`Intercept(tfine)`$coef
  beta <- coefs$smterms$`ff(X)`$coef
  residuals <- model$residuals
  
  return(list("alpha" = alpha, "beta" = beta, "residuals" = residuals))
  
}

A = Sys.time()

for (b in 1:1000) {
  


print(paste0("Simulation ", b, " of 1000"))


# read in data
  x.fd = eval.fd(tfine, Data[[3]][[b]]$x.fd)
  y.fd = eval.fd(tfine, Data[[3]][[b]]$yplusE)
  
pffr.data = data.frame("ID" = 1:100)
pffr.data$Y = t(y.fd)
pffr.data$X = t(x.fd)

## fit the model on sim data

fit = pffr(Y ~ ff(X, limits = "s<=t", xind = tfine,
                  splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                    k = c(20, 20))), yind = tfine, data = pffr.data, bs.int = list(bs = 'ps', k = 20, m = c(2,1)))


Fit.List[[b]] = get.info(model = fit)



}


# save sim data to save time 

details = "pffr fit, k = c(20, 20), m = c(2, 1)"


sim_results <- list("Details" = details, "Results" = Fit.List)
saveRDS(sim_results, "pffr Simulation Results.rds")

A - Sys.time()

