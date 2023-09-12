library(rtdists)
library(pmwg)
library(R.matlab)
library(dplyr)
rm(list=ls())
data <- readMat("LBA_Forstmann_scanner.mat")
data <- data$data


tmp <- data.frame(matrix(ncol = 5, nrow = 0))
names <- c("response", "rt", "cond","scan", "subject")
colnames(tmp) <- names

for (i in 1:20){
  x<-as.data.frame(cbind(unlist(data[[1]][i]),unlist(data[[2]][i]),unlist(data[[3]][i]),unlist(data[[4]][i])))
  x$subject <- rep(i,nrow(x))
  names(x)<-names
  
  tmp<-bind_rows(tmp,x)
  
}
rownames(tmp)<-NULL
data <- tmp
data <- data[data$rt>.15,]


lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  if (any(data$rt < x["t0"])) {
    return(-1e10)
  }
  if (sample){
    data$rt=NA
    data$resp = NA
  }
  
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$cond]
  
  if (sample) {
    out <- rtdists::rLBA(n = nrow(data),
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = x[c("v1", "v2")],
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    data$rt <- out$rt
    data$response <- out$resp
    
  } else {
    out <- rtdists::dLBA(rt = data$rt,
                         response = data$response,
                         A = x["A"],
                         b = bs,
                         t0 = x["t0"],
                         mean_v = list(x["v1"],x[ "v2"]),
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    bad <- (out < 1e-10) | (!is.finite(out))
    out[bad] <- 1e-10
    out <- sum(log(out))
  }
  if (sample){return(data)}
  if (!sample){return(out)}
}

ll=function(x,data,sample=FALSE){
  ### here i pull in the parameters for each component of the joint model and rename them for the LL. For safety, you could also do this by using grepl() and only taking those with certain values (i.e. "in.b1" or "out.b1"). 
  xIn<-x[1:7]
  xOut<-x[8:14]
  names(xIn)<-names(xOut)<- c("A","b1","b2","b3","v1","v2","t0")
  
  if (!sample){
    out.In = lba_loglike(xIn,data[data$scan==1,],sample=FALSE)
    out.Out = lba_loglike(xOut,data[data$scan==0,],sample=FALSE)
    
    ### check if either returned bad values (from bad parameters) and if so, return a bad outcome, otherwise, sum.
    if(out.In > -1e+10 && out.Out > -1e+10){
      #### this is the linking part. We can sum the two outcomes
      out <- sum(out.In,out.Out)
    }    else {
      out <- -1e10
    }
    return(out)
  }
  
  if (sample){
    ### as the data is structured similarly and our likleihood is setup to return the same structure with generated rt and resp, this part is simple. Other LLs may require further manipulation
    In <- lba_loglike(xIn,data[data$scan==1,],sample=TRUE)
    Out <- lba_loglike(xOut,data[data$scan==0,],sample=TRUE)
    data <- rbind(In, Out)
    return(data)
  }
}

pars <- c("in.b1", "in.b2", "in.b3", "in.A", "in.v1", "in.v2", "in.t0",
          "out.b1", "out.b2", "out.b3", "out.A", "out.v1", "out.v2", "out.t0")

priors <- list(
  theta_mu_mean = rep(0, length(pars)),
  theta_mu_var = diag(rep(1, length(pars)))) 


sampler <- pmwgs(
  data = data,
  pars = pars,
  prior = priors,
  ll_func = ll
)

sampler <- init(sampler)

sampled <- run_stage(sampler, stage = "burn",iter = 100, particles = 100, n_cores = 16)


sampled <- run_stage(sampled, stage = "adapt",iter = 10000, particles = 100, n_cores =16)


sampled <- run_stage(sampled, stage = "sample", iter = 1000, particles = 100, n_cores = 16)

save.image("bev.RData")



