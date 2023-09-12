setwd ("")	# Change to your working directory


########################## READ IN THE DATA

load ("mentalRotation.RData")
source ("lba-math.r")

nPpn = length (unique (datErp[,11]))
nCon = 5
nStim = 2

Ppn = unique (datErp[,11])
datErp = datErp[datErp[,13]<7000,]
datErp[datErp[,12]>100,12] = datErp[datErp[,12]>100,12]-100

NeuroData = list ()
for (i in 1:nPpn)
{
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

tmp=do.call(rbind,NeuroData)

########################## IMPORTANT FUNCTIONS

# log density of subject-level parameters under group-level truncated
# normal distributions.
log.dens.hyper = function (x, hyper)
{
	out = 0  
	for (p in names(x))
	{
		out = out + dtnorm (x[p], hyper[paste (p, "mu", sep = ".")], hyper[paste (p, "sigma", sep = ".")], 0, Inf, log = TRUE)
	}
	out
}

# log density for the data (RT & response) under subject-level parameters.
# Assumes data$S and data$R have identical levels !!!
log.dens.like = function (x, data, par.names, robust = robust, posdrift = posdrift)
{
	names(x) = par.names
  
	# get one parameter per accumulator 
	nresps = 2
	As = rep (x["A"], nresps)
	bs = c(x["bS"] + x["A"], x["bM"] + x["A"])
	ss = rep (1, nresps)
	offset = x["offset"] - 10
	
	# Loop over conditions, stimuli and responses to get likelihoods
	out = numeric (length (data$RT))
	for (i in 1:nCon)
	{
		vSC = x[paste("vSC", (i-1)*45, sep = ".")]	# 5 same rotation correct drifts
		vSE = x[paste("vSE", (i-1)*45, sep = ".")]	# 5 same rotation error drifts
		vMC = x[paste("vMC", (i-1)*45, sep = ".")]	# 5 mirror rotation correct drifts
		vME = x[paste("vME", (i-1)*45, sep = ".")]	# 5 mirror rotation error drifts

		for (j in 1:nStim)
		{
			vs = switch (j, c(vSC, vME), c(vSE, vMC))

			for (k in 1:nresps)
			{
				n1 = c(k, -k+3) # index of 1st rate parameter
				isSR =  data$C==i & data$S==j & data$R==k
				out[isSR] = log (pmax (n1PDFfixedt0 (data$RT[isSR] - x["t0"], x0max = As[n1], chi = bs[n1], drift = vs[n1], s = ss[n1], posdrift = posdrift, robust = robust), 1e-10)) + 
				dnorm (data$ERP23[isSR], mean = vs[k] * x["scale23"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP34[isSR], mean = vs[k] * x["scale34"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP45[isSR], mean = vs[k] * x["scale45"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP56[isSR], mean = vs[k] * x["scale56"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP67[isSR], mean = vs[k] * x["scale67"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP78[isSR], mean = vs[k] * x["scale78"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP89[isSR], mean = vs[k] * x["scale89"] + offset, sd = x["sd"], log = TRUE) +
				dnorm (data$ERP90[isSR], mean = vs[k] * x["scale90"] + offset, sd = x["sd"], log = TRUE)
			}
		}
	}
	sum (out)
}

# log density of the subject-level parameters under the group-level
# truncated normal distribution, plus log density of the parameters
# for the group-level distribution under the priors.
log.dens.hyper.and.prior = function (theta, phi, prior)
{
	sum ((dtnorm (theta, phi[1], phi[2], 0, Inf, log = TRUE))) +
		(dtnorm (phi[1], prior$mu[1], prior$mu[2], 0, Inf, log = TRUE)) +
		(dgamma (phi[2], prior$sigma[1], prior$sigma[2], log = TRUE))
}

# Runs the MCMC step using DE proposals, for parameters defined at the
# individual subject level. "i" is the chain to update
crossover = function (i, pars, use.theta, use.like, data, par.names, hyper, robust, posdrift, rp)
{
	# Combine old likelihood and old prior to get old posterior density
	use.weight = use.like[i] + log.dens.hyper (use.theta[i,], hyper[i,])
	gamma = 2.38 / sqrt (2 * length (pars)) # step size
	index = sample (c(1:nChains)[-i], 2, replace = F) # pick two other chains
	theta = use.theta[i,]	# old parameters					
	theta[pars] = use.theta[i,pars] + # DE step
		gamma * (use.theta[index[1], pars] - use.theta[index[2], pars]) + runif (1, -rp, rp)
	weight = log.dens.hyper (theta, hyper[i,]) # likelihood of proposal params under group-level.
	# Only bother looking up likelihood of data under new proposals, if those proposal
	# parameters have non-zero likelihood under the group-level.
	if (is.finite (weight))
	{
		like = log.dens.like (theta, data, par.names = par.names, robust = robust, posdrift = posdrift)
	}
	else
	{
		like = 0
	}
	weight = like + weight # New posterior, without prior (which hasn't changed).
	if (is.na (weight)) weight = -Inf
	if (runif (1) < exp (weight - use.weight))
	{					
		use.theta[i,] = theta
		use.like[i] = like
	}
	c(use.like[i], use.theta[i,])
}

# Runs the MCMC step using DE proposals, for group-level parameters.
crossover_hyper = function (i, pars, use.theta, use.phi, prior, rp)
{
	use.weight = log.dens.hyper.and.prior (use.theta[i,], use.phi[i,pars], prior)
	gamma = 2.38 / sqrt (2 * length (pars))
	index = sample (c(1:nChains)[-i], 2, replace = F)
	phi = use.phi[i,]
	phi[pars] = use.phi[i,pars] + gamma * (use.phi[index[1], pars] - use.phi[index[2], pars]) + runif (1, -rp, rp)
	weight = suppressWarnings (log.dens.hyper.and.prior (use.theta[i,], phi[pars], prior))
	if (is.na (weight)) weight = -Inf
	if (runif (1) < exp (weight - use.weight)) use.phi[i,] = phi
	use.phi[i,]
}


########################## INITIALISE MODEL

begin = date()

posdrift = TRUE	# Use truncated normal (strictly positive) drift rate distributions?
robust = FALSE	# Use the robustness checks in lba-math.r functions? Safer, but slower.

nChains = 32	# usually 2-3 x number of subject level parameters or more.
vSC = paste ("vSC", seq (0, 180, 45), sep = ".")
vSE = paste ("vSE", seq (0, 180, 45), sep = ".")
vMC = paste ("vMC", seq (0, 180, 45), sep = ".")
vME = paste ("vME", seq (0, 180, 45), sep = ".")
theta.names = c("A", "bS", "bM", vSC, vSE, vMC, vME, "scale23", "scale34", "scale45", "scale56", "scale67", "scale78", "scale89", "scale90", "offset", "sd", "t0") # actually b here is B = b-A
sensible.start.points = c(2, rep (1, 2), rep (3, 5), rep (1, 5), rep (3, 5), rep (1, 5), rep (1, 8), 8, 5, 0.2)		# Mean start points for chains
names (sensible.start.points) = theta.names
nPars = length (theta.names)
nHpars = 2 * nPars
nmc = 2000
rp = .001		# This is the random perturbation in the DE-MCMC proposals.

theta = array (NA, c(nChains, nPars, nPpn, nmc))	# parameters
phi = array (NA, c(nChains, nHpars, nmc))			# group-level parameters
weight = array (-Inf, c(nmc, nChains, nPpn))		# likelihood of parameter vector 
colnames (theta) = theta.names
colnames (phi) = paste (rep (theta.names, each = 2), c("mu", "sigma"), sep = ".")

# Truncated normal priors for location (mu) parameters and gamma priors for stdev (sigma) parameters. 
prior = list()
prior$A = prior$bS = prior$bM = list (mu = c(2, 1), sigma = c(1, 1))
for (i in 1:length (vSC)) {prior[[vSC[i]]] = list (mu = c(3, 1), sigma = c(1, 1))}
for (i in 1:length (vSE)) {prior[[vSE[i]]] = list (mu = c(1, 1), sigma = c(1, 1))}
for (i in 1:length (vMC)) {prior[[vMC[i]]] = list (mu = c(3, 1), sigma = c(1, 1))}
for (i in 1:length (vME)) {prior[[vME[i]]] = list (mu = c(1, 1), sigma = c(1, 1))}
prior$scale23 = prior$scale34 = prior$scale45 = prior$scale56 = prior$scale67 = prior$scale78 = prior$scale89 = prior$scale90 = list (mu = c(1, 1), sigma = c(1, 1))
prior$offset = list (mu = c(8, 2), sigma = c(1, 1))
prior$sd = list (mu = c(5, 1), sigma = c(1, 1))
prior$t0 = list (mu = c(0.2, .1), sigma = c(1, 3))

# Generate some fairly widely-spread random start points.
for (i in 1:nChains)
{
	# First get start points for subject-level parameters.
	for (j in 1:nPpn)
	{
		while (weight[1,i,j]==-Inf)
		{ # Make sure the parameters are valid.
			theta[i,,j,1] = rtnorm (n = nPars, mean = sensible.start.points, sd = sensible.start.points / 10, 0, Inf)
			weight[1,i,j] = log.dens.like (theta[i,,j,1], data = NeuroData[[j]], par.names = theta.names, robust = robust, posdrift = posdrift)
		}
	}
	# Now for group-level parameters.
	for (p in theta.names)
	{
		phi[i,paste (p, "mu", sep = "."), 1] = rtnorm (n = 1, mean = sensible.start.points[p], sd = sensible.start.points[p] / 10, lower = 0)
		phi[i,paste (p, "sigma", sep = "."), 1] = rtnorm (n = 1, mean = sensible.start.points[p] / 10, sd = sensible.start.points[p] / 20, lower = 0)
	}   
}


########################## FIT MODEL

# Run the sampling.
for (i in 2:nmc)
{
	cat ("\n ", i, "  ")

	# First sample new group-level parameters.
	phi[,,i]=phi[,,i-1]
	for (p in theta.names)
	{
		which.theta = match (x = p, table = theta.names)
		which.phi = match (x = paste (p, c("mu", "sigma"), sep = "."), table = colnames (phi))
		phi[,,i] = t(sapply (1:nChains, crossover_hyper, pars = which.phi, 
			use.theta = theta[,which.theta,,i-1], use.phi = phi[,,i], prior = prior[[p]], rp = rp))
	}

	# Now sample new subject-level parameters.
	for(j in 1:nPpn)
	{
		cat(j)
		temp = t(sapply (1:nChains, crossover, pars = 1:nPars, use.theta = theta[,,j,i-1], use.like = weight[i-1,,j], 
			data = NeuroData[[j]], hyper = phi[,,i], par.names = theta.names, robust = robust, posdrift = posdrift, rp = rp))
		weight[i,,j] = temp[,1]
		theta[,,j,i] = temp[,2:(nPars + 1)]
	}
}

Means = matrix (, nPpn, nPars); for (i in 1:nPpn) {for (j in 1:nPars) {Means[i,j] = mean (theta[,j,i,])}}; colnames (Means) = theta.names
PhiMeans = c(); for (j in 1:nHpars) {PhiMeans[j] = mean (phi[,j,])}; names (PhiMeans) = colnames (phi)


########################## SYNTHETIC DATA

# A function to generate fake data that match the likelihood
# function, to be used for testing, and for generating the
# posterior predictive data. For one subject!
synthesise.data = function (x, n, posdrift = TRUE)
{
	# x should be a vector of parameters in same format as used
	# by log.dens.like(), and n should inform sample size.
	data = list (S = NULL, C = NULL, R = NULL, RT = NULL, ERP23 = NULL, ERP34 = NULL, ERP45 = NULL, ERP56 = NULL, ERP67 = NULL, ERP78 = NULL, ERP89 = NULL, ERP90 = NULL)
    bs = c(x["bS"] + x["A"], x["bM"] + x["A"])
    
	# Loop over the 5 rotation conditions
	for (i in 1:nCon)
	{
		vSC = x[paste("vSC", (i-1)*45, sep = ".")] # 5 same rotation correct drifts
		vSE = x[paste("vSE", (i-1)*45, sep = ".")] # 5 same rotation error drifts
		vMC = x[paste("vMC", (i-1)*45, sep = ".")] # 5 mirror rotation correct drifts
		vME = x[paste("vME", (i-1)*45, sep = ".")] # 5 mirror rotation error drifts

		# Loop over same and mirror stimuli
		for (j in 1:nStim)
		{
			vs = switch (j, c(vSC, vME), c(vSE, vMC))
			tmp = rlba (n = n[j,i], b = bs, A = x["A"], vs = vs, s = 1, t0 = x["t0"], posdrift = posdrift)
			tmp23 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale23"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp34 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale34"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp45 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale45"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp56 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale56"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp67 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale67"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp78 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale78"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp89 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale89"] + (x["offset"] - 10), sd = x["sd"]) 
			tmp90 = rnorm (n = n[j,i], mean = vs[tmp$resp] * x["scale90"] + (x["offset"] - 10), sd = x["sd"]) 
			data$S = c(data$S, rep (j, length (tmp$rt)))
			data$C = c(data$C, rep (i, length (tmp$rt)))
			data$R = c(data$R, tmp$resp)
			data$RT = c(data$RT, tmp$rt)
			data$ERP23 = c(data$ERP23, tmp23)
			data$ERP34 = c(data$ERP34, tmp34)
			data$ERP45 = c(data$ERP45, tmp45)
			data$ERP56 = c(data$ERP56, tmp56)
			data$ERP67 = c(data$ERP67, tmp67)
			data$ERP78 = c(data$ERP78, tmp78)
			data$ERP89 = c(data$ERP89, tmp89)
			data$ERP90 = c(data$ERP90, tmp90)
		}
	}
	data$S = factor (data$S)
	data$R = factor (data$R)
	data$C = factor (data$C)
	data.frame (data)
}

# Generate some posterior predictive data to assess model fit.
burnin = nmc - 1000		# Guess for burn-in level.

n.posterior = 25 # Number of parameter samples from posterior distribution.
pp.data = list ()
for (s in 1:nPpn)
{
	cat(s, " ")
	ns = table (NeuroData[[s]]$S, NeuroData[[s]]$C)
	pp.data[[s]] = list (S = NULL, C = NULL, R = NULL, RT = NULL, ERP23 = NULL, ERP34 = NULL, ERP45 = NULL, ERP56 = NULL, ERP67 = NULL, ERP78 = NULL, ERP89 = NULL, ERP90 = NULL)
	for (i in round (seq (from = burnin, to = nmc, length.out = n.posterior)))
	{
		tmp = synthesise.data (x = theta[sample (nChains, 1),,s,i], n = ns, posdrift = posdrift)
		pp.data[[s]]$S = c(pp.data[[s]]$S, tmp$S)
		pp.data[[s]]$C = c(pp.data[[s]]$C, as.character (tmp$C))
		pp.data[[s]]$R = c(pp.data[[s]]$R, tmp$R)
		pp.data[[s]]$RT = c(pp.data[[s]]$RT, tmp$RT)
		pp.data[[s]]$ERP23 = c(pp.data[[s]]$ERP23, tmp$ERP23)
		pp.data[[s]]$ERP34 = c(pp.data[[s]]$ERP34, tmp$ERP34)
		pp.data[[s]]$ERP45 = c(pp.data[[s]]$ERP45, tmp$ERP45)
		pp.data[[s]]$ERP56 = c(pp.data[[s]]$ERP56, tmp$ERP56)
		pp.data[[s]]$ERP67 = c(pp.data[[s]]$ERP67, tmp$ERP67)
		pp.data[[s]]$ERP78 = c(pp.data[[s]]$ERP78, tmp$ERP78)
		pp.data[[s]]$ERP89 = c(pp.data[[s]]$ERP89, tmp$ERP89)
		pp.data[[s]]$ERP90 = c(pp.data[[s]]$ERP90, tmp$ERP90)
	}
	pp.data[[s]] = data.frame (pp.data[[s]])
}
names (pp.data) = names (NeuroData)

end = date(); begin; end
save.image ("Fit.RData")


########################## CONVERGENCE PLOTS

jpeg (file = "Convergence.jpg", width = 1500, height = 1500)

par (mfcol = c(12, 12), mar = c(2, 2, 1, 1))	# Adjust mfcol to have prod more than 2*nHpars

for (i in 1:nHpars)
{
	plot (y = as.vector (phi[,i,1:nmc]), x = rep (1:nmc, each = nChains), type = "p", pch = 16, cex = .5, xlab = "", ylab = "")  
	mtext (side = 3, dimnames (phi)[[2]][i])
	hist (phi[,i,burnin:nmc], main = "", xlab = "", ylab = "")
}

dev.off ()


########################## QUANTILES

# RTs
qps = seq (.1, .9, .2)	# Change this for finer or coarser granularity
tmp = lapply (NeuroData, function (x) tapply (x$RT, list (x$S==x$R, x$S, x$C), quantile, prob = qps))
for (i in 1:nPpn) {for (j in 1:2) {for (k in 1:nStim) {for (l in 1:nCon) {if (length (tmp[[i]][j,k,l][[1]]) == 0) tmp[[i]][j,k,l][[1]] = rep (NA, length (qps))}}}}
q = array (unlist (tmp), dim = c(length (qps), 2, nStim, nCon, length (tmp)), dimnames = list (qps, c("Error", "Correct"), dimnames (tmp[[1]])[[2]], dimnames (tmp[[1]])[[3]], names (tmp)))
tmp = lapply (pp.data, function (x) tapply (x$RT, list (x$S==x$R, x$S, x$C), quantile, prob = qps))
for (i in 1:nPpn) {for (j in 1:2) {for (k in 1:nStim) {for (l in 1:nCon) {if (length (tmp[[i]][j,k,l][[1]]) == 0) tmp[[i]][j,k,l][[1]] = rep (NA, length (qps))}}}}
pp.q = array (unlist (tmp), dim = c(length (qps), 2, nStim, nCon, length (tmp)), dimnames = list (qps, c("Error", "Correct"), dimnames (tmp[[1]])[[2]], dimnames (tmp[[1]])[[3]], names (tmp)))

# Pcs
tmp = lapply (NeuroData, function (x) tapply (x$S==x$R, list (x$S, x$C), mean))
if (length (unlist (tmp))!=(nStim * nCon * length (tmp))) stop ("Some missing data in accuracy calculations for real data.")
p = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)), dimnames = list (names (tmp[[1]]), names (tmp[[2]]), names (tmp)))
tmp = lapply (pp.data, function (x) tapply (x$S==x$R, list (x$S, x$C), mean))
if (length (unlist (tmp))!=(nStim * nCon * length (tmp))) stop ("Some missing data in accuracy calculations for pp data.")
pp.p = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)), dimnames = list (names (tmp[[1]]), names (tmp[[2]]), names (tmp)))

# ERPs
ERPs = c("ERP23", "ERP34", "ERP45", "ERP56", "ERP67", "ERP78", "ERP89", "ERP90")
mERP = pp.mERP = array (NA, dim = c(8, 2, 5, 9), dimnames = list (2:9, c("Same", "Mirror"), seq (0, 180, 45), 1:9))

for (i in 1:8)
{
	tmp = lapply (NeuroData, function (x) tapply (unlist (x[ERPs[i]]), list (x$S, x$C), mean))
	mERP[i,,,] = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)))
	tmp = lapply (pp.data, function (x) tapply (unlist (x[ERPs[i]]), list (x$S, x$C), mean))
	pp.mERP[i,,,] = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)))
}


########################## POSTERIOR PREDICTIVE PLOTS

pdf ("PP.pdf")

m = matrix (1:20, nrow = 5)
m = rbind (0, cbind (0, m[,1], 0, m[,2], 0, m[,3], 0, m[,4], 0), 0)
layout (mat = m, widths = c(.6, 1, .6, 1, .6, 1, .6, 1, .1), heights = c(.1, rep (1, nCon), .7))
par (mar = rep (0, 4), cex.axis = 0.8)

# First plot average data, dots are data, lines are fit.
av.q = apply (q, 1:4, mean, na.rm = T); av.pp.q = apply (pp.q, 1:4, mean, na.rm = T)
av.p = apply (p, 1:2, mean, na.rm = T); av.pp.p = apply (pp.p, 1:2, mean, na.rm = T)
av.mERP = apply (mERP, 1:3, mean, na.rm = T); av.pp.mERP = apply (pp.mERP, 1:3, mean, na.rm = T)

for (j in 1:nStim)
{
	for (k in 1:nCon)
	{
		matplot (x = av.q[,,j,k], y = cbind (qps * (1 - av.p[j,k]), qps * av.p[j,k]), type = "p",
			pch = 16, col = c("red", "green"), cex = .5, axes = F, xlab = "", ylab = "", ylim = c(0, 1), xlim = c(0.3, 2.5))
		axis (side = 1, labels = k==nCon); axis (side = 2, labels = NA)
		mtext (side = 2, line = 0.6, (k-1)*45)
		if (k==3) mtext (side = 2, line = 1.8, dimnames(av.mERP)[[2]][j])		
		if (k==nCon) mtext (side = 1, line = 4, c("                        RT (sec.)", "Group                        ")[j])
		matlines (x = av.pp.q[,,j,k], y = cbind (qps * (1 - av.pp.p[j,k]), qps * av.pp.p[j,k]), type = "l", lty = 1, col = c("red", "green"))
	}
}

# Now plot average ERP data, dots are data, lines are fit.
for (j in 1:nStim)
{
	for (k in 1:nCon)
	{
		plot (x = (1:8)+0.5, y = av.mERP[,j,k], type = "p", pch = 16, col = c("green"), cex = .5, axes = F, xlab = "", ylab = "", xlim = c(1, 9), ylim = c(0, 5))
		lines (x = (1:8)+0.5, y = av.pp.mERP[,j,k], type = "l", lty = 1, col = c("green"))
		axis (side = 1, at = 1:9, labels = if (k==5) seq (200, 1000, 100) else NA); axis (side = 2, labels = NA)
		mtext (side = 2, line = 0.6, (k-1)*45)
		if (k==3) mtext (side = 2, line = 1.8, dimnames(av.mERP)[[2]][j])		
		if (k==nCon) mtext (side = 1, line = 4, c("                      ERP epoch (ms.)", "Group             ")[j])
	}
}

dev.off ()


########################## SUMMARY EFFECT PLOTS

ci = function (x) {if (sd (x) == 0) {0} else {t.test(x)$conf.int[1:2]}}
tmp = lapply (NeuroData, function (x) tapply (x$RT, list (x$S, x$C), median))
MRT = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)), dimnames = list (dimnames (tmp[[1]])[[1]], dimnames (tmp[[1]])[[2]], names (tmp)))
tmp = lapply (pp.data, function (x) tapply (x$RT, list (x$S, x$C), median))
pp.MRT = array (unlist (tmp), dim = c(nStim, nCon, length (tmp)), dimnames = list (dimnames (tmp[[1]])[[1]], dimnames (tmp[[1]])[[2]], names (tmp)))
av.MRT = apply (MRT, 1:2, mean); av.pp.MRT = apply (pp.MRT, 1:2, mean)
ci.MRT = apply (MRT, 1:2, ci)
ci.p = apply (p, 1:2, ci)
ci.ERP = apply (mERP, 1:3, ci); mci.ERP = apply (ci.ERP, c(1, 3, 4), mean)

pdf ("Plots/Summary Effects4.pdf")

layout (matrix (1:4, 2, 2))
par (mar = c (4, 4, 1, 1), cex.axis = 0.7)

plot (x = NA, axes = F, xlab = "", ylab = "", xlim = c(1, 10.5), ylim = c(0, 2))
axis (1, at = c(1:5, seq (6.5, 10.5, 1)), labels = rep (seq (0, 180, 45), 2)); axis (2)
for (i in 1:nCon)
{
	lines (x = i, y = av.MRT[1,i], type = "p", pch = 16, cex = 1.5, col = i)
	lines (x = i+5.5, y = av.MRT[2,i], type = "p", pch = 16, cex = 1.5, col = i)
	arrows (x0 = i, y0 = ci.MRT[1,1,i], y1 = ci.MRT[2,1,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
	arrows (x0 = i+5.5, y0 = ci.MRT[1,2,i], y1 = ci.MRT[2,2,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
}
mtext (side = 1, line = 2, "Same                       Mirror")
mtext (side = 1, line = 3, "Condition")
mtext (side = 2, line = 2, "Median RT")
lines (x = 1:5, y = av.pp.MRT[1,], lty = 2)
lines (x = seq (6.5, 10.5, 1), y = av.pp.MRT[2,], lty = 2)
legend (x = 5, y = 2, pch = c(16, NA), lty = c(NA, 2), legend = c("Data", "Model"))

plot (x = NA, axes = F, xlab = "", ylab = "", xlim = c(1, 10.5), ylim = c(0.5, 1))
axis (1, at = c(1:5, seq (6.5, 10.5, 1)), labels = rep (seq (0, 180, 45), 2)); axis (2)
for (i in 1:nCon)
{
	lines (x = i, y = av.p[1,i], type = "p", pch = 16, cex = 1.5, col = i)
	lines (x = i+5.5, y = av.p[2,i], type = "p", pch = 16, cex = 1.5, col = i)
	arrows (x0 = i, y0 = ci.p[1,1,i], y1 = ci.p[2,1,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
	arrows (x0 = i+5.5, y0 = ci.p[1,2,i], y1 = ci.p[2,2,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
}
mtext (side = 1, line = 2, "Same                       Mirror")
mtext (side = 1, line = 3, "Condition")
mtext (side = 2, line = 2, "Mean Pc")
lines (x = 1:5, y = av.pp.p[1,], lty = 2)
lines (x = seq (6.5, 10.5, 1), y = av.pp.p[2,], lty = 2)

plot (x = NA, axes = F, xlab = "", ylab = "", xlim = c(1, 10), ylim = c(0, 5))
mtext (side = 1, line = 2, "ms.")
mtext (side = 2, line = 2, "ERP Wave (Same)")
axis (side = 1, at = 1:10, labels = c(seq (200, 1000, 100), "Av")); axis (side = 2)
legend (x = 6, y = 5, fill = 1:5, legend = seq (0, 180, 45))
for (i in 1:nCon)
{
	lines (x = 1:8+.5, y = av.mERP[,1,i], type = "p", pch = 16, cex = 1.5, col = i)
	lines (x = 1:8+.5, y = av.pp.mERP[,1,i], type = "l", cex = 1.5, col = i)
	arrows (x0 = 9+.2*i-.1, y0 = mci.ERP[1,1,i], y1 = mci.ERP[2,1,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
	# for (l in 1:8)
	# {
		# arrows (x0 = l+.2*i-.1, y0 = ci.ERP[1,l,1,i], y1 = ci.ERP[2,l,1,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
	# }		
}

Link = matrix (, 3, 8)
for (i in 1:8) {Link[,i] = quantile (phi[,i*2+45,], prob = seq (.25, .75, .25))}
plot (x = NA, axes = F, xlab = "", ylab = "", xlim = c(1, 10), ylim = c(0, 5))
mtext (side = 1, line = 2, "ms.")
mtext (side = 2, line = 2, "ERP Wave (Mirror)")
axis (side = 1, at = 1:10, labels = c(seq (200, 1000, 100), "Av")); axis (side = 2)
text (x = 3, y = 4.2, pos = 4, labels = "Linking parameter")
for (i in 1:nCon)
{
	lines (x = 1:8+.5, y = av.mERP[,2,i], type = "p", pch = 16, cex = 1.5, col = i)
	lines (x = 1:8+.5, y = av.pp.mERP[,2,i], type = "l", cex = 1.5, col = i)
	arrows (x0 = 9+.2*i-.1, y0 = mci.ERP[1,2,i], y1 = mci.ERP[2,2,i], code = 3, angle = 90, length = 0.05, lwd = 1.5, col = i)
	for (l in 1:8)
	{
		lines (x = rep (l+.2, 2), y = 3+c(0, Link[2,l]), type = "l")
		lines (x = rep (l+.8, 2), y = 3+c(0, Link[2,l]), type = "l")
		lines (x = c(l+.2, l+.8), y = rep (3, 2), type = "l")
		lines (x = c(l+.2, l+.8), y = rep (3+Link[2,l], 2), type = "l")
		arrows (x0 = l+.5, y0 = 3+Link[1,l], y1 = 3+Link[3,l], code = 3, angle = 90, length = 0.05, lwd = 1.5)
	}		
}

dev.off ()

