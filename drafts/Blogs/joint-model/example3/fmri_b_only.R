library(rtdists)
library(fmri)
library(pmwg)
library(tibble)
rm(list=ls())


generateOnsets <- function(scans=1500, TR=1.5, trialDuration=5, jitterOptions=c(0,1,2,3), addNullTrials=TRUE) {
  times <- seq(1, scans*TR, TR)
  trialDuration <- 5  # in seconds
  jitterOptions <- c(0, 1, 2, 3)
  trialOnsets <- seq(0, max(times), trialDuration)
  trialOnsets <- trialOnsets + sample(jitterOptions, length(trialOnsets), replace=TRUE)
  
  if(addNullTrials) {
    # 10% null trials, in an actual experiment you'd pseudo-randomize the null trial onsets to ensure they're not too predictable
    # but in this simulation it doesn't really matter so just remove every 10th onset
    trialOnsets <- trialOnsets[-seq(10, length(trialOnsets), 10)]
  }
  
  # half is easy, half is hard
  easyOnsets <- sort(sample(trialOnsets, size=length(trialOnsets)/2, replace=FALSE))
  hardOnsets <- setdiff(trialOnsets, easyOnsets)
  
  return(list('easy'=easyOnsets, 'hard'=hardOnsets))
}

generateDM <- function(onsets, TR, scans) {
  # convolve onsets with HRF - use 'fmri' package here, I'm too lazy
  easy <- fmri.stimulus(scans=scans, onsets=onsets[['easy']]/TR, TR=TR)
  hard <- fmri.stimulus(scans=scans, onsets=onsets[['hard']]/TR, TR=TR)
  
  # design matrix - first column is intercept
  dm <- cbind(rep(1, length(easy)), easy, hard)
  
  return(dm)
}

# function to 'simulate' fmri data in single ROI
simulatefMRI <- function(beta0, beta1, beta2, dm, noisesd=0.02, plot=FALSE) {
  # The simulated data is a weighted sum of both event types, plus some noise (let's simplify and assume Gaussian noise [this is not the case in real data])
  # so we need parameters weighting the influence of both event types (easy and hard)
  simulatedData <- as.vector(c(beta0, beta1, beta2) %*% t(dm)) + rnorm(nrow(dm), 0, noisesd)
  
  if(plot) {
    par(mfrow=c(2,1))
    plot(dm[,2], type='l', col=1, xlab='Scan N', ylab='Expected signal (a.u.)', xlim=c(1,100))
    lines(dm[,3], col=2)
    legend('topleft', legend=c('Easy', 'Hard'), col=1:2, lty=c(1,1), bty='n')
    plot(simulatedData, type='l', col=1, xlab='Scan N', ylab='Expected signal (a.u.)', xlim=c(1,100), main='simulated signal')
    title(sub=paste(c('beta0', 'beta1', 'beta2'), round(c(beta0, beta1, beta2),2), sep = '=', collapse='  '))
  }
  return(simulatedData)
}

simulateBehavior <- function(b, vMean, vEasy, vHard, t0, onsets) {
  # use onsets to get number of easy and hard trials
  
  # easy <- rdiffusion(length(onsets$easy), v=vEasy, a=a, t0=t0)
  # hard <- rdiffusion(length(onsets$hard), v=vHard, a=a, t0=t0)
  easyDrift<-list(vMean-vEasy, vMean+vEasy)
  hardDrift<-list(vMean-vHard, vMean+vHard)
  
  easy <- rLBA(n = length(onsets$easy), mean_v = easyDrift, A = 1, b=b, t0=t0, sd_v = c(1,1))
  hard <- rLBA(n = length(onsets$hard), mean_v = hardDrift, A = 1, b=b, t0=t0, sd_v = c(1,1))
  
  easy$onset <- onsets$easy
  hard$onset <- onsets$hard
  easy$trialType <- 'easy'
  hard$trialType <- 'hard'
  
  simulatedData <- rbind(easy, hard)
  return(simulatedData[order(simulatedData$onset),])
}


nSubs <- 32

# DDM parameters
b <- runif(n=nSubs, 1, 1)
vMean <- runif(n=nSubs, 2, 4)
vEasy <- rnorm(n=nSubs, 2, 1)
vHard <- vEasy / runif(n=nSubs, 1.5, 2)
t0 <- runif(n=nSubs, 0.1, 0.3)

# fMRI parameters
beta0 <- runif(n=nSubs, 10, 30)  # intercept, not so interesting
# let's assume our BOLD responses are linearly dependent on drift rates, such that beta = w*v + some noise
w <- runif(n=nSubs, 0.3, 0.6)
beta1 <- vEasy*w + rnorm(nSubs, 0.05, 0.05)
beta2 <- vHard*w + rnorm(nSubs, 0.05, 0.05)
noisesd <- runif(nSubs, 0.02, 0.1)  # let's also estimate noise


# For each subject, generate a set of onsets, make the design matrix (dm), and then simulate both the fMRI and the behavioral data

fMRIdata=matrix(nrow=1500, ncol = nSubs)
behavioralData <- NULL

for (i in 1:nSubs){
  onsets <- generateOnsets(scans=1500, TR=1.5)
  dm <- generateDM(onsets, scans=1500, TR=1.5)
  fMRIdata[,i] <- simulatefMRI(beta0[i], beta1[i], beta2[i], noisesd[i], dm=dm, plot=FALSE)
  Data <- simulateBehavior(b[i], vMean[i], vEasy[i], vHard[i], t0[i], onsets)
  Data$subject <- i
  behavioralData <- rbind(behavioralData, Data) 
}

DM<-as.data.frame(dm)
data<-NULL
for (i in 1:nSubs){
  df <- tibble(
    subject = i,
    behavioralData = list(
      tibble(behavioralData[behavioralData$subject==i,])
    ),
    fmriData = list(
      tibble(fMRIdata[,i])
    ),
    dm = list(
      tibble(DM)
    )
  )
  data<-rbind(data,df)
}


llBehavior <- function(A, b, vMean, vEasy, vHard, t0, behavioralData) {
  if (any(behavioralData$rt < t0 )) {
    return(-1e10)
  }
  idxEasy <- behavioralData$trialType=='easy'
  idxHard <- behavioralData$trialType=='hard'
  
  easyDrift<-list(vMean-vEasy, vMean+vEasy)
  hardDrift<-list(vMean-vHard, vMean+vHard)
  
  behavioralData$like[idxEasy] <- dLBA(rt=behavioralData$rt[idxEasy],resp=behavioralData$response[idxEasy], mean_v = easyDrift, A = A, b=b, t0=t0, sd_v = c(1,1), silent = TRUE)
  behavioralData$like[idxHard] <- dLBA(rt=behavioralData$rt[idxHard],resp=behavioralData$response[idxHard], mean_v = hardDrift, A = A, b=b, t0=t0, sd_v = c(1,1), silent = TRUE)
  
 
  out <- sum(log(pmax(behavioralData$like, 1e-10)))
  
  if (is.numeric(out)){
    return(out)
  }else{
    out <- -1e+10
    return(out)
  }
}

llfMRI <- function(beta0, beta1, beta2, noisesd, fMRIdata, dm) {
  
  expectedSignal <- as.vector(c(beta0, beta1, beta2) %*% t(dm))
  likes <- dnorm(unlist(fMRIdata)-expectedSignal, sd=noisesd)
  
  
  out <- sum(log(pmax(likes, 1e-10)))
  
  if (is.numeric(out)){
    return(out)
  }else{
    out <- -1e+10
    return(out)
  }
}

ll <- function(x,data,sample=FALSE){
  x<- exp(x)
  out <- llBehavior(A=x["A"],b=x["b"], vMean=x["vMean"], vEasy = x["vEasy"], vHard = x["vHard"], t0 = x["t0"], behavioralData=as.data.frame(data$behavioralData))
  #out.neural <- llfMRI(beta0 = x["beta0"], beta1 = x["beta1"], beta2 = x["beta2"], noisesd = x["noise"], fMRIdata = as.data.frame(data$fmriData), dm = as.data.frame(data$dm))
  
  
  
    return(out)
}




pars <- c("A","b","vMean", "vEasy", "vHard", "t0")

priors <- list(
  theta_mu_mean = c(0,0,1,1,1,-2),
  theta_mu_var = diag(c(1,1,1,1,1,0.5))) 


sampler <- pmwgs(
  data = data,
  pars = pars,
  prior = priors,
  ll_func = ll
)

sampler <- init(sampler)
#save.image("fmri.RData")

sampled <- run_stage(sampler, stage = "burn",iter = 100, particles = 50, n_cores = 16)
save.image("fmri.RData")


sampled <- run_stage(sampled, stage = "adapt",iter = 5000, particles = 100, n_cores = 16)
save.image("fmri.RData")


sampled <- run_stage(sampled, stage = "sample", iter = 1000, particles = 100, n_cores = 16)

save.image("fmri.RData")



