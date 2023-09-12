##### fMRI and Behavioural data
library(rtdists)
library(fmri)
library(lme4)

### task
#### 3 conditions - easy, medium and hard
#### 3 blocks per condition, 20 trials per block (3x3x20 = 180 Trials)
#### 2 choices - 1 & 2
#### b.data = sub, cond, stim, resp, rt
#### n.data = sub, cond, event, bold


#### onsets = all times of trials, durations = length of each
#### pars = a1,a2,b1,b2,cc

#### 15 seconds per trial
conds = c("easy", "medium", "hard")
n.cond = length(conds)
n.trials = 20
n.blocks = n.conds * 3
n.subj=10
names = c("subject","resp","rt")
data = data.frame(matrix(NA, ncol = length(names), nrow = (n.trials*n.subj*n.cond))) #empty data frame
names(data)=names
data$condition = rep(1:n.cond,times = n.trials) #filling in condition
data$subject = rep(1:n.subj, each = n.trials*n.cond) #filling in subjects


lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]
  
  if (sample) {
    data$rt=NA
    data$resp = NA
    for (i in 1:nrow(data)){
      out <- rtdists::rLBA(n = 1,
                           A = x["A"],
                           b = bs[i],
                           t0 = x["t0"],
                           mean_v = x[c("v1", "v2")],
                           sd_v = c(1, 1),
                           distribution = "norm",
                           silent = TRUE)
      data$rt[i] <- out$rt
      data$resp[i] <- out$resp
    }
    return(data)
  } else {
    if (any(data$rt < x["t0"])) {
      return(-1e10)
    }
    out <- rtdists::dLBA(rt = data$rt,
                         response = data$resp,
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
    return(out)
  }
}


parameter.names=c("b1","b2","b3", "A","v1","v2","t0",
                  "beta.0","beta.1","beta.2","beta.3", "within.sd")
n.parameters=length(parameter.names)

ptm <- array(dim = n.parameters, dimnames = list(parameter.names)) #an empty array where i will put parameter values
pts2 <- array(dim=c(n.parameters,n.parameters), dimnames = list(parameter.names,parameter.names)) #an empty array where i will put the covariance matrix

ptm[1:n.parameters]=c(0.1,0.3,0.5,0.4,1.2,0.3,-2,
                      3,4,3,2,0.5) 
vars = abs(ptm)/10 #off diagonal correlations are done as absolute/10
sigmaC = matrix(.1,nrow=length(ptm),ncol=length(ptm)) ###std dev correlation on diagonal - you might think this should be corr = 1, but it's actually the standard deviation 
diag(sigmaC)=sqrt(vars)
sigmaC[8:12,8:12]<-0.4
sigmaC[1:3,1:3]<-0.6
sigmaC[1:3,7]<-0.4
sigmaC[7,1:3]<-0.4
sigmaC[8:11,1:3]<-0.5
sigmaC[1:3,8:11]<-0.5


sigmaC <- sdcor2cov(sigmaC)

subj_random_effects <- t(mvtnorm::rmvnorm(n.subj,mean=ptm,sigma=sigmaC))

for (i in 1:n.subj){
  tmp<- lba_loglike(subj_random_effects[,i],sample=TRUE,data=data[data$subject==i,])
  data$rt[data$subject==i]=tmp$rt
  data$resp[data$subject==i]=tmp$resp
}
head(data)

#T=n.trials * 3 * n.blocks

# trial <- c(0,1,0)
# block <- rep(trial, 20)


## in this example, fMRI only collected for one block of data
emptytrial <- rep(0,10)
tasktrial <- rep(1,10)
EVstim <- rep(c(emptytrial,tasktrial,emptytrial),20)
T=length(EVstim)*3


tmp <- which(EVstim>0)
baseonsets <- seq(from =min(tmp), to = max(tmp)-9, by = 30 )
EV1 <- fmri.stimulus(scans=T, onsets=baseonsets,durations=rep(10, n.trials), type="canonical", par =
                       c(6,12,0.9,0.9,0.15))

EV2 <- fmri.stimulus(scans=T, onsets=baseonsets+600,durations=rep(10, n.trials), type="canonical", par =
                      c(4,6,0.7,0.9,0.15))

EV3 <- fmri.stimulus(scans=T, onsets=baseonsets+1200,durations=rep(10, n.trials), type="canonical", par =
                      c(2,2,0.7,0.9,0.10))


plot(1:T, EV1,type="n", ylab="BOLD Response", xlab="Volume Number")
lines(1:T,EV1)
lines(1:T,EV2)
lines(1:T,EV3)


generateSubjectVoxelData <- function(i) {
  Y <- ptm[8] + EV1*ptm[9] + EV2*ptm[10] + EV3*ptm[11] +
    subj_random_effects[8,i] + EV1*subj_random_effects[9,i] + EV2*subj_random_effects[10,i] + EV3*subj_random_effects[11,i] + EV3*beta3.betweensub[i] +
    rnorm(T, 0,pnorm(subj_random_effects[12,i]))                                                                                     
  return(Y)
}

# subj_random_effects[,i]
# beta0.fixed <- 3
# beta1.fixed <- 2
# beta2.fixed <- 5
# beta3.fixed <- 4
# within.sd <- 0.5

allsubjs <- lapply(1:n, generateSubjectVoxelData)
Y.g <- unlist(allsubjs)
subject <- unlist(lapply(1:n,rep,T))
subject <- as.factor(subject)
EV1.vec <- rep(EV1, n)
EV2.vec <- rep(EV2, n)
EV3.vec <- rep(EV3, n)
voxeldat <- data.frame(subject, Y.g, EV1.vec, EV2.vec, EV3.vec)

voxeldat$source <- "fMRI"
data$source <- "behaviour"


tmp <- data %>%
  split(.$subject) %>%
  tibble(subject = names(.), behaviour = .) 

tmp2 <- voxeldat %>%
  split(.$subject) %>%
  tibble(subject = names(.), fMRI = .) 

data <- merge(tmp,tmp2, by = "subject")
data$subject<- as.numeric(data$subject)


tmp <- data[data$subject==1,]
tmp=do.call(rbind,tmp$fMRI)



lba_loglike <- function(x, data, sample = FALSE) {
  x <- exp(x)
  bs <- x["A"] + x[c("b1", "b2", "b3")][data$condition]
  
  if (sample) {
    data$rt=NA
    data$resp = NA
    for (i in 1:nrow(data)){
    out <- rtdists::rLBA(n = 1,
                         A = x["A"],
                         b = bs[i],
                         t0 = x["t0"],
                         mean_v = x[c("v1", "v2")],
                         sd_v = c(1, 1),
                         distribution = "norm",
                         silent = TRUE)
    data$rt[i] <- out$rt
    data$resp[i] <- out$resp
    }
    return(data)
  } else {
    if (any(data$rt < x["t0"])) {
      return(-1e10)
    }
    out <- rtdists::dLBA(rt = data$rt,
                         response = data$resp,
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
    return(out)
  }
}

fMRI_like = function(){
  if(sample){
    
  }
  else {
    
  }
}

LL = function(){
  data.fMRI <- do.call(rbind,data$fMRI)
  data.behaviour <- do.call(rbind,data$behaviour)

  if (!sample){
    out.fMRI = lba_loglike(xIn,data[data$scan==1,],sample=FALSE)
    out.behaviour = lba_loglike(xOut,data[data$scan==0,],sample=FALSE)
    
    ### check if either returned bad values (from bad parameters) and if so, return a bad outcome, otherwise, sum.
    if(out.fMRI > -1e+10 && out.behavour > -1e+10){
      #### this is the linking part. We can sum the two outcomes
      out <- sum(out.fMRI,out.behaviour)
    }    else {
      out <- -1e10
    }
    return(out)
  }
  
  if (sample){
    ### as the data is structured similarly and our likleihood is setup to return the same structure with generated rt and resp, this part is simple. Other LLs may require further manipulation
    fMRI.data <- lba_loglike(xIn,data[data$scan==1,],sample=TRUE)
    behaviour.data <- lba_loglike(xOut,data[data$scan==0,],sample=TRUE)
    
    return(data)
  } 
}





