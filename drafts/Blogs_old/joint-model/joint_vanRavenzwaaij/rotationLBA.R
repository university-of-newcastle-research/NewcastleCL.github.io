rm(list=ls())

library(dplyr)
library(pmwg)
library(rtdists)

##load in the data using van Ravenzwaaij code
load ("mentalRotation.RData")

nPpn = length (unique (datErp[,11]))
nCon = 5
nStim = 2

Ppn = unique (datErp[,11])
datErp = datErp[datErp[,13]<7000,]
datErp[datErp[,12]>100,12] = datErp[datErp[,12]>100,12]-100

NeuroData = list ()
for (i in 1:nPpn){
  NeuroData[[i]] = list (RT = NULL, C = NULL, S = NULL, R = NULL, ERP = NULL, ERP23 = NULL, ERP34 = NULL, ERP45 = NULL, ERP56 = NULL, ERP67 = NULL, ERP78 = NULL, ERP89 = NULL, ERP90 = NULL)
  NeuroData[[i]]$RT = datErp[datErp[,11]==Ppn[i],13]/1000
  NeuroData[[i]]$C = (datErp[datErp[,11]==Ppn[i],12]-1)%%5+1										# 1 = 0, 2 = 45, 3 = 90, 4 = 135, 5 = 180
  NeuroData[[i]]$S = (datErp[datErp[,11]==Ppn[i],12]<21)+1										# 1 = same, 2 = mirror
  Correct = datErp[datErp[,11]==Ppn[i],14]														# 1 = correct, 0 = incorrect
  NeuroData[[i]]$R = (NeuroData[[i]]$S * Correct) + ((-NeuroData[[i]]$S + 3) * (-Correct + 1))	# 1 = same, 0 = mirror
  NeuroData[[i]]$ERP = apply (datErp[datErp[,11]==Ppn[i],3:10], 1, mean)							# ERP mean from 200 to 1000 ms
  for (j in 3:10) {NeuroData[[i]][[j+3]] = datErp[datErp[,11]==Ppn[i],j]}	# Individual ERP bins
  # NeuroData[[i]]$S = factor (NeuroData[[i]]$S)
  # NeuroData[[i]]$R = factor (NeuroData[[i]]$R)
  # NeuroData[[i]]$C = factor (NeuroData[[i]]$C)
  #NeuroData[[i]]$Subject = rep(i,length(NeuroData[[i]]$RT))
  NeuroData[[i]] = data.frame(NeuroData[[i]])
}

for (i in 1:nPpn){
  NeuroData[[i]]$subject <- rep(i,length(NeuroData[[i]]$RT))
}
data=do.call(rbind,NeuroData)
data <- data[data$RT < 3,]

ll=function(x,data,sample=FALSE){
  #i only exp 1-4 here as the two coefficients could be pos or neg
  x[1:4]=exp(x[1:4])
  if (any(data$rt < c(x["t0"]) )) {
    return(-1e10)
  }
  v1=v2=out=numeric(nrow(data))
  for (c in unique(data$C)){
    for (s in unique(data$S)){
      use = data$C == c & data$S == s
      v1[use]<- ifelse(s==1, x["v.mean"] + x["ERP89"]*data$ERP89[use], 
                       x["v.mean"] - x["ERP89"]*data$ERP89[use])
      v2[use]<- ifelse(s==2, x["v.mean"] + x["ERP89"]*data$ERP89[use], 
                       x["v.mean"] - x["ERP89"]*data$ERP89[use] )
    }
  }
  if(any(v1<0)){
    return(-1e10)
  }
  if(any(v2<0)){
    return(-1e10)
  }
  if (sample){
    data$rt <- NA
    data$response <- NA
    for (i in 1:nrow(data)){
      tmp=rLBA(n=1,
               A=A,b=bs,t0=t0,mean_v=list(v1[i],v2[i]),sd_v=list(1,1),
               distribution = "norm", silent = TRUE)
      data$RT[i] <- tmp$rt
      data$R[i] <- tmp$resp
    }
    return(data)
  }else{
    out=dLBA(rt=data$RT,response = data$R,
             A=x["A"],b=x["b"],t0=x["t0"],mean_v=list(v1,v2),sd_v=list(1,1),
             distribution = "norm", silent = TRUE)
    
  }
  return(sum(log(pmax(out, 1e-10))))
}

pars <- c("A", "b","t0", "v.mean","ERP89")      

priors <- list(
  theta_mu_mean = rep(0, length(pars)),
  theta_mu_var = diag(rep(1, length(pars)))) 

sampler <- pmwgs(
  data = data,
  pars = pars,
  prior = priors,
  ll_func = ll
)

start_points <- list(
  mu = c(log(c(1,.5,0.2,3)),0.2),
  sig2 = diag(rep(.1, length(pars)))
)



sampler <- init(sampler, start_mu = start_points$mu,
                start_sig = start_points$sig2)


sampled <- run_stage(sampler, stage = "burn",iter = 1000, particles = 100, n_cores = 16)

sampled <- run_stage(sampled, stage = "adapt", iter = 5000,particles = 100, n_cores = 16)

sampled <- run_stage(sampled, stage = "sample", iter = 1000,particles = 50, n_cores = 16)


save.image("rotation_covariate.Rdata")
