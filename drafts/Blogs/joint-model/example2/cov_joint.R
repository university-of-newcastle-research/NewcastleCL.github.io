library(rtdists)
library(pmwg)
library(dplyr)
rm(list=ls())

load ("mentalRotation.RData")

#the code below comes from https://osf.io/2r7bv/
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
  NeuroData[[i]]$S = factor (NeuroData[[i]]$S)
  NeuroData[[i]]$R = factor (NeuroData[[i]]$R)
  NeuroData[[i]]$C = factor (NeuroData[[i]]$C)
  NeuroData[[i]] = data.frame (NeuroData[[i]])
}
for (i in 1:nPpn){
  NeuroData[[i]]$subject <- rep(i,length(NeuroData[[i]]$RT))
}
tmp=do.call(rbind,NeuroData)


ERPavg <- tmp %>% group_by(subject, C, S) %>% summarise(erp89 = mean(ERP89),
                                                        erp90 = mean(ERP90)) %>% ungroup()

data<-left_join(tmp[,c(1:4,14)],ERPavg, by = c("subject", "C", "S"))
data$C<-as.character(data$C)
data$R<-as.character(data$R)
data$S<-as.character(data$S)
data <- data[data$RT<4,]

ll=function(x,data,sample=FALSE){
  #i only exp 1-4 here as the two coefficients could be pos or neg
  x[1:4]=exp(x[1:4])
  if (any(data$RT < c(x["t0"]) )) {
    return(-1e10)
  }
  A=x["A"]
  t0=x["t0"]
  bs = x["b"]
  v1=v2=out=numeric(nrow(data))
  for (c in unique(data$C)){
    for (s in unique(data$S)){
      use = data$C == c & data$S == s
      v1[use]<- ifelse(s==1, x["v.mean"] + x["ERP89"]*data$erp89[use] + x["ERP90"]*data$erp90[use], 
                       x["v.mean"] - x["ERP89"]*data$erp89[use] - x["ERP90"]*data$erp90[use])
      v2[use]<- ifelse(s==2, x["v.mean"] + x["ERP89"]*data$erp89[use] + x["ERP90"]*data$erp90[use], 
                       x["v.mean"] - x["ERP89"]*data$erp89[use] - x["ERP90"]*data$erp90[use])
    }
  }
  # if(any(v1) < 0){
  #   return(-1e10)
  # }
  # if(any(v2) < 0){
  #   return(-1e10)
  # }
  if (sample){
    data$RT <- NA
    data$R <- NA
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
             A=A,b=bs,t0=t0,mean_v=list(v1,v2),sd_v=list(1,1),
             distribution = "norm", silent = TRUE)
    
  }
  return(sum(log(pmax(out, 1e-10))))
}

pars <- c("A", "b","t0", "v.mean","ERP89","ERP90")      


priors <- list(
  theta_mu_mean = c(0,0,-2,1,-1,-1),
  theta_mu_var = diag(rep(1, length(pars)))) 


sampler <- pmwgs(
  data = data,
  pars = pars,
  prior = priors,
  ll_func = ll
)

sampler <- init(sampler)

sampled <- run_stage(sampler, stage = "burn",iter = 200, particles = 100, n_cores = 16)


sampled <- run_stage(sampled, stage = "adapt",iter = 10000, particles = 100, n_cores =16)
save.image("cov.RData")


sampled <- run_stage(sampled, stage = "sample", iter = 1000, particles = 100, n_cores = 16)

save.image("cov.RData")



