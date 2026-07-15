library(fda)
library(tidyverse)
library(reshape)
library(VectorCodeR)
library(gsignal)

setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")

basis = create.bspline.basis(c(0,1), 20, 4)

hip.fd = fd(gait[,,1], basis)
knee.fd = fd(gait[,,2], basis)


t = seq(0,1,0.01)
e.hip = eval.fd(t, hip.fd) 
e.knee = eval.fd(t, knee.fd)

# vector coding

vector_coding = unwrap(phase_angle(e.hip, e.knee))


VC.Mean = rowMeans(vector_coding)
VC.SD = apply(vector_coding, 1, sd)
t2 = seq(0,1, length = 100)

P1 = ggplot()+
  geom_ribbon(aes(x = t2, ymin = VC.Mean-VC.SD, ymax = VC.Mean+VC.SD), fill = "red", alpha = 0.4)+
  geom_line(aes(x= t2, y = VC.Mean), lty = 1)+
  theme_light()+
  labs(x = "t", y = "Coupling Angle")
P1

ggsave("HFLMHM - Figure 7.jpg", dpi = 300, width = 10, height = 10, units = 'cm')

# continuous relative phase

# normalise signals

e.hippad = rbind(e.hip[81:101,],e.hip, e.hip[1:20,])
e.kneepad = rbind(e.knee[81:101,],e.knee, e.knee[1:20,])


hip = apply(e.hippad, 2, function(X){
  
  X - min(X) - (max(X)-min(X))/2
  
})

knee = apply(e.kneepad, 2, function(X){
  
  X - min(X) - (max(X)-min(X))/2
  
})

# hilbert transform
h.hip = apply(hip, 2, function(X){Im(hilbert(X))})
h.knee = apply(knee, 2, function(X){Im(hilbert(X))})

# calculate CRP
crp.wrapped = atan2(h.hip*knee-hip*h.knee, hip*knee+h.hip*h.knee) # corrected

CRP = unwrap(crp.wrapped)*(180/pi)

par(mfrow = c(1,1))
matplot(CRP, type = 'l', xlab = 't', ylab = "Continous Relative Phase", main = "Corrected")

# get mean and SD (variability)

CRP.Mean = rowMeans(CRP)
CRP.SD = apply(CRP, 1, sd)

P2 = ggplot()+
  geom_ribbon(aes(x = t, ymin = CRP.Mean[22:122]-CRP.SD[22:122], ymax = CRP.Mean[22:122]+CRP.SD[22:122]), fill = "blue", alpha = 0.4)+
  geom_line(aes(x= t, y = CRP.Mean[22:122]), lty = 1)+
  theme_light()+
  labs(x = "t", y = "Continuous Relative Phase")
P2

ggsave("HFLM - Figure 8.jpg", dpi = 300, width = 10, height = 10, units = 'cm')

