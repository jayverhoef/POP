library(POP)
library(Matrix)
library(viridis)
library(classInt)
library(rgeos)
library(batchmeans)

data(POPhexagons)
data(AKpoly)

# extract the data from POPhexagons
	ptDF = POPhexagons@data

# these data objects were created from POPhexagons using the preliminaries.R script
# they are needed for using MCMC to fit the models
data(nb_list)
data(W)
data(Dmat)
data(rhogrid)
data(detsig)
data(SigiList)


#-------------------------------------------------------------------------------
# Prepare for MCMC
#-------------------------------------------------------------------------------

# extract the observed data vector, and which are sampled and missing
y = ptDF$shipdays 
# find which cells have no effort, nor any of their neighbors
ind_no_nei_effort = sapply(nb_list,function(x){all(is.na(y[x]))}) & is.na(y)
# set these to zero effort
y[ind_no_nei_effort] = 0
# create indicator of those with certain effort or zero effort
ind_samp = which(!is.na(y))
# the rest are missing - unknown effort
ind_miss = which(is.na(y))

#check creation of hard zeros
# get the centroids of the polygons for faster plotting
trueCentroids = gCentroid(POPhexagons, byid = TRUE)
# create breaks
fixedBreaks = c(0:2, 10, 30, 70, max(y, na.rm = TRUE) + 1 + 1e-10) - .5
# get colors for classes
ffb = classIntervals(y, style = 'fixed', fixedBreaks = fixedBreaks)
pal = viridis(length(fixedBreaks) - 1)
pal[1] = 'maroon3'
ffbColours = findColours(ffb, pal)
# plot the data
plot(AKpoly[1:100,], col = 'grey50')
plot(trueCentroids, col = ffbColours, pch = 19, add = TRUE, cex = .6)
plot(AKpoly[1:100,], add = TRUE, border = 'black')
text(650000,1850000,'Ship-days', cex = 2.7, col = 'darkorange1')
text(650000,1600000,'Zeros Added', cex = 2.7, col = 'darkorange1')


#-------------------------------------------------------------------------------
# Do MCMC
#-------------------------------------------------------------------------------

# bounds for truncating the multivariate normal distribution
trunc_bound = 6

# define some alternative values for the data distribution in the hierarchical model
dpoi_mu = function(x, mu, sigma, log = TRUE)
{
	dpois(x, lambda = mu, log = log)
}

dgamma_mu_sigma = function(x, mu, sigma, log = TRUE)
{
	dgamma(x, shape = mu^2/sigma^2, scale = sigma^2/mu, log = log)
}

dnbinom_mu_sigma = function(x, mu, sigma, log = TRUE)
{
	dnbinom(x, mu = mu, size = sigma, log = log)
}

#set a random number seed so results are reproducible
set.seed(1001)
#first run MCMC for fixed sigma, sigma_CAR, and rho to let z and beta mix
start.time = Sys.time()
MCMCout = MCMC(response = y, datadist = dpoi_mu, 
  beta0ini = .3, sigmaini = .5, sigCARini = 1.2, 
	rhoini = rhogrid[66], zini = rep(0, times = length(y)),
	rhogrid = rhogrid, detsig = detsig, SigiList = SigiList, 
	ind_samp = ind_samp, ind_miss = ind_miss, 
	sample_sigma = FALSE, sample_sigCAR = FALSE,
	sample_rho = FALSE, use_trunc = TRUE, trunc_bound = trunc_bound,
	nMCMC = 1000, n_z_iters = 100, n_thin = 20, samp_sd = .01, miss_sd = .01,
	beta0_sd = .1, sigma_sd = .01, sigCAR_sd = .02)
end.time = Sys.time()
# check the time of the run
difftime(end.time, start.time, units = 'mins')

# check acceptance rates
# for sampled hexagons (including introduced zeros)
sum(MCMCout$accept_samp)/length(MCMCout$accept_samp)
# for unsampled hexagons
sum(MCMCout$accept_miss)/length(MCMCout$accept_miss)
# for beta0
sum(MCMCout$accept_beta0)/length(MCMCout$accept_beta0)


# get the number of MCMC samples from previous output
# so MCMC can start up again from the final set of values
n = length(MCMCout$store_beta0)
#set a random number seed so results are reproducible
set.seed(1002)
# treat this run as a burn-in
start.time = Sys.time()
MCMCout = MCMC(response = y, datadist = dpoi_mu, 
	beta0ini = MCMCout$store_beta0[n], sigmaini = MCMCout$store_sigma[n],
	sigCARini = MCMCout$store_sigCAR[n], 
	rhoini = MCMCout$store_rho[n],
	zini = MCMCout$store_zmat[,n], rhogrid = rhogrid, 
	detsig = detsig, SigiList = SigiList, 
	ind_samp = ind_samp, ind_miss = ind_miss, sample_sigma = FALSE, 
	sample_sigCAR = TRUE, sample_rho = TRUE, 
	use_trunc = TRUE, trunc_bound = trunc_bound,
  nMCMC = 1000, n_z_iters = 50, n_thin = 50, samp_sd = .02, miss_sd = .02,
	beta0_sd = .1, sigma_sd = 1, sigCAR_sd = .02, rho_indxpm = 3)
end.time = Sys.time()
# check the time of the run
difftime(end.time, start.time, units = 'mins')

# check acceptance rates
# for sampled hexagons (including introduced zeros)
sum(MCMCout$accept_samp)/length(MCMCout$accept_samp)
# for unsampled hexagons
sum(MCMCout$accept_miss)/length(MCMCout$accept_miss)
# for beta0
sum(MCMCout$accept_beta0)/length(MCMCout$accept_beta0)
# for sigma (data model variance)
sum(MCMCout$accept_sigma)/length(MCMCout$accept_sigma)
# for sigma_CAR (CAR error variance)
sum(MCMCout$accept_sigCAR)/length(MCMCout$accept_sigCAR)
# for rho (CAR autocorrelation parameter)
sum(MCMCout$accept_rho)/length(MCMCout$accept_rho)

#check the limits of z
max(MCMCout$store_zmat)
min(MCMCout$store_zmat)
#plot some traces of sampled values
plot(MCMCout$store_sigCAR, type = 'l')
plot(log(MCMCout$store_sigma), type = 'l')
plot(MCMCout$store_rho, type = 'l')
plot(MCMCout$store_beta0, type = 'l')
plot(MCMCout$store_zmat[14,], type = 'l')
plot(MCMCout$store_zmat[5001,], type = 'l')

# again, so we can start sampling with final values from burn-in
n = length(MCMCout$store_beta0)
#set a random number seed so results are reproducible
set.seed(1003)
# final run
start.time = Sys.time()
MCMCout = MCMC(response = y, datadist = dpoi_mu, 
	beta0ini = MCMCout$store_beta0[n], sigmaini = MCMCout$store_sigma[n],
	sigCARini = MCMCout$store_sigCAR[n], 
	rhoini = MCMCout$store_rho[n],
	zini = MCMCout$store_zmat[,n], rhogrid = rhogrid, 
	detsig = detsig, SigiList = SigiList, 
	ind_samp = ind_samp, ind_miss = ind_miss, sample_sigma = FALSE, 
	sample_sigCAR = TRUE, sample_rho = TRUE, 
	use_trunc = TRUE, trunc_bound = trunc_bound,
  nMCMC = 1000, n_z_iters = 50, n_thin = 50, samp_sd = .02, miss_sd = .03,
	beta0_sd = .05, sigma_sd = 3, sigCAR_sd = .02, rho_indxpm = 3)
end.time = Sys.time()
difftime(end.time, start.time, units = 'mins')

# check acceptance rates
sum(MCMCout$accept_samp)/length(MCMCout$accept_samp)
sum(MCMCout$accept_miss)/length(MCMCout$accept_miss)
sum(MCMCout$accept_beta0)/length(MCMCout$accept_beta0)
sum(MCMCout$accept_sigCAR)/length(MCMCout$accept_sigCAR)
sum(MCMCout$accept_rho)/length(MCMCout$accept_rho)

#check the limits of z
max(MCMCout$store_zmat)
min(MCMCout$store_zmat)
#plot some traces of sampled values
plot(MCMCout$store_sigCAR, type = 'l')
plot(MCMCout$store_rho, type = 'l')
plot(MCMCout$store_beta0, type = 'l')
plot(MCMCout$store_zmat[14,], type = 'l')
plot(MCMCout$store_zmat[5001,], type = 'l')
	
# Convergence Diagnostics
logfit = 0*MCMCout$store_zmat
for(i in 1:length(MCMCout$store_beta0)) logfit[,i] = MCMCout$store_beta0[i] + MCMCout$store_zmat[,i]
ESS_all = apply(logfit,1,ess)
min(ESS_all)
BM_all = bmmat(t(logfit))
median(BM_all[,2])
max(BM_all[,2])

mode = function(x, bw = "nrd0")
{
	densout = density(x, bw = bw)
	densout$x[densout$y == max(densout$y)]	
}

# get the fitted means per MCMC sample and then get averages and 
# standard deviations per hexagon across all MCMC samples
temp = matrix(NA, nrow = dim(MCMCout$store_zmat)[1], ncol = dim(MCMCout$store_zmat)[2])
for(i in 1:dim(MCMCout$store_zmat)[2])
	temp[,i] = exp(MCMCout$store_beta0[i] + MCMCout$store_zmat[,i])
mode_post = apply(temp,1,mode)
sd_post = apply(temp,1,sd)
zsd_post = apply(MCMCout$store_zmat,1,sd)

#check the result graphically
fixedBreaks = c(0, 1, 2, 10, 30, 70, max(mode_post, na.rm = TRUE) + 1e-10) - 1e-12
f10 = classIntervals(mode_post, style = 'fixed', fixedBreaks = fixedBreaks)
pal = viridis(length(f10$brks) -1)
f10Colours = findColours(f10, pal)
plot(AKpoly[1:100,], col = 'grey20')
plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, pch = 19, cex = .7)
plot(AKpoly[1:100,], add = TRUE, border = 'black')
text(650000,1950000,'Ship-Day', cex = 2.2, col = 'white')
text(650000,1710000,'Posterior', cex = 2.2, col = 'white')
text(650000,1450000,' ', cex = 2.2, col = 'white')
addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
    colors = pal, cex = 1.1, printFormat = "6.0")

# create objects from list items
MCMC_eff_accept_samp = MCMCout$accept_samp
MCMC_eff_accept_miss = MCMCout$accept_miss
MCMC_eff_accept_beta0 = MCMCout$accept_beta0
MCMC_eff_accept_sigCAR = MCMCout$accept_sigCAR
MCMC_eff_accept_rho = MCMCout$accept_rho
MCMC_eff_zmat = MCMCout$store_zmat
MCMC_eff_beta0 = MCMCout$store_beta0
MCMC_eff_sigCAR = MCMCout$store_sigCAR
MCMC_eff_rho = MCMCout$store_rho

#store the results in the package so they can be quickly accessed for graphics and 
#further processing
setwd(paste0('/home/jay/data/2019_papers/POP',
  '/POP_package/POP/data'))
save(MCMC_eff_zmat, file = 'MCMC_eff_zmat.rda')
save(MCMC_eff_beta0, file = 'MCMC_eff_beta0.rda')
save(MCMC_eff_sigCAR, file = 'MCMC_eff_sigCAR.rda')
save(MCMC_eff_rho, file = 'MCMC_eff_rho.rda')
save(MCMC_eff_accept_samp, file = 'MCMC_eff_accept_samp.rda')
save(MCMC_eff_accept_miss, file = 'MCMC_eff_accept_miss.rda')
save(MCMC_eff_accept_beta0, file = 'MCMC_eff_accept_beta0.rda')
save(MCMC_eff_accept_sigCAR, file = 'MCMC_eff_accept_sigCAR.rda')
save(MCMC_eff_accept_rho, file = 'MCMC_eff_accept_rho.rda')

