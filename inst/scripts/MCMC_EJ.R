library(POP)
library(Matrix)
library(viridis)
library(classInt)
library(rgeos)
library(batchmeans)
library(nabor)

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
# for northern fur seals
y = ptDF$EJ 

#get MCMC results from fitting effort surface
data(MCMC_eff_zmat)
data(MCMC_eff_beta0)

#get posterior distribution of effort per cell, and mean per cell
MCMC_effdays = matrix(NA, nrow = dim(MCMC_eff_zmat)[1], 
	ncol = dim(MCMC_eff_zmat)[2])
for(i in 1:dim(MCMC_eff_zmat)[2])
	MCMC_effdays[,i] = exp(MCMC_eff_zmat[,i] + MCMC_eff_beta0[i])
MCMC_eff_mean = apply(MCMC_effdays,1,mean)
# create an indicator of at least one day of effort from posterior mean of ship-days
ind_eff = !is.na(ptDF$shipdays)  |  MCMC_eff_mean >= 1 

# if sufficient effort, but species not observed, set grid cell to zero
y[ind_eff & is.na(y)] = 0
# create indicator of those with data, including newly added zeros
ind_samp = which(!is.na(y))
# the rest are missing - unknown effort
ind_miss = which(is.na(y))


#check creation of hard zeros
# get the centroids of the polygons for faster plotting
trueCentroids = gCentroid(POPhexagons, byid = TRUE)
# create breaks
fixedBreaks = c(0:2, 6, 12, 25, max(y, na.rm = TRUE) + 1 + 1e-10) - .5
# get colors for classes
ffb = classIntervals(y, style = 'fixed', fixedBreaks = fixedBreaks)
pal = viridis(length(fixedBreaks) - 1)
pal[1] = 'maroon3'
ffbColours = findColours(ffb, pal)
# plot the data
plot(AKpoly[1:100,], col = 'grey50')
plot(trueCentroids, col = ffbColours, pch = 19, add = TRUE, cex = .6)
plot(AKpoly[1:100,], add = TRUE, border = 'black')
text(650000,1850000,'Steller Sea Lions', cex = 2.7, col = 'darkorange1')
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
set.seed(2001)

# Initial run to get values populated
start.time = Sys.time()
#first run MCMC for fixed sigma, sigma_CAR, and rho to let z and beta mix
MCMCout = MCMC(response = y, datadist = dnbinom_mu_sigma, 
  beta0ini = -2.5, sigmaini = 1, sigCARini = 1.5, 
	rhoini = rhogrid[66], zini = rep(0, times = length(y)),
	rhogrid = rhogrid, detsig = detsig, SigiList = SigiList,  eff_mat_offset = MCMC_effdays,
	ind_samp = ind_samp, ind_miss = ind_miss, 
	sample_sigma = FALSE, sample_sigCAR = FALSE,
	sample_rho = FALSE, use_trunc = TRUE, trunc_bound = trunc_bound,
	nMCMC = 1000, n_z_iters = 100, n_thin = 20, samp_sd = .01, miss_sd = .02,
	beta0_sd = .1, sigma_sd = .1, sigCAR_sd = .02)
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
set.seed(2002)
# treat this run as a burn-in
start.time = Sys.time()
MCMCout = MCMC(response = y, datadist = dnbinom_mu_sigma, 
	beta0ini = MCMCout$store_beta0[n], sigmaini = MCMCout$store_sigma[n],
	sigCARini = MCMCout$store_sigCAR[n], 
	rhoini = MCMCout$store_rho[n],
	zini = MCMCout$store_zmat[,n], rhogrid = rhogrid, 
	detsig = detsig, SigiList = SigiList, eff_mat_offset = MCMC_effdays,
	ind_samp = ind_samp, ind_miss = ind_miss, sample_sigma = TRUE, 
	sample_sigCAR = TRUE, sample_rho = TRUE, 
	use_trunc = TRUE, trunc_bound = trunc_bound,
  nMCMC = 1000, n_z_iters = 50, n_thin = 50, samp_sd = .02, miss_sd = .03,
	beta0_sd = .1, sigma_sd = .2, sigCAR_sd = .02, rho_indxpm = 5)
end.time = Sys.time()
difftime(end.time, start.time, units = 'mins')
sum(MCMCout$accept_samp)/length(MCMCout$accept_samp)
sum(MCMCout$accept_miss)/length(MCMCout$accept_miss)
sum(MCMCout$accept_beta0)/length(MCMCout$accept_beta0)
sum(MCMCout$accept_sigma)/length(MCMCout$accept_sigma)
sum(MCMCout$accept_sigCAR)/length(MCMCout$accept_sigCAR)
sum(MCMCout$accept_rho)/length(MCMCout$accept_rho)

# again, so we can start sampling with final values from burn-in
n = length(MCMCout$store_beta0)
#set a random number seed so results are reproducible
set.seed(2003)
# final run
start.time = Sys.time()
MCMCout = MCMC(response = y, datadist = dnbinom_mu_sigma, 
	beta0ini = MCMCout$store_beta0[n], sigmaini = MCMCout$store_sigma[n],
	sigCARini = MCMCout$store_sigCAR[n], 
	rhoini = MCMCout$store_rho[n],
	zini = MCMCout$store_zmat[,n], rhogrid = rhogrid, 
	detsig = detsig, SigiList = SigiList, eff_mat_offset = MCMC_eff_mean,
	ind_samp = ind_samp, ind_miss = ind_miss, sample_sigma = TRUE, 
	sample_sigCAR = TRUE, sample_rho = TRUE, 
	use_trunc = TRUE, trunc_bound = trunc_bound,
  nMCMC = 1000, n_z_iters = 50, n_thin = 50, samp_sd = .05, miss_sd = .07,
	beta0_sd = .1, sigma_sd = .2, sigCAR_sd = .02, rho_indxpm = 5)
end.time = Sys.time()
difftime(end.time, start.time, units = 'mins')
sum(MCMCout$accept_samp)/length(MCMCout$accept_samp)
sum(MCMCout$accept_miss)/length(MCMCout$accept_miss)
sum(MCMCout$accept_beta0)/length(MCMCout$accept_beta0)
sum(MCMCout$accept_sigma)/length(MCMCout$accept_sigma)
sum(MCMCout$accept_sigCAR)/length(MCMCout$accept_sigCAR)
sum(MCMCout$accept_rho)/length(MCMCout$accept_rho)

#check the limits of z
max(MCMCout$store_zmat)
min(MCMCout$store_zmat)
#plot some traces of sampled values
plot(MCMCout$store_sigma, type = 'l')
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
fixedBreaks = c(0, .1, .2, .5, 1, 5, max(mode_post, na.rm = TRUE) + 1e-10) - 1e-12
f10 = classIntervals(mode_post, style = 'fixed', fixedBreaks = fixedBreaks)
pal = viridis(length(f10$brks) -1)
f10Colours = findColours(f10, pal)
plot(AKpoly[1:100,], col = 'grey20')
plot(trueCentroids, col = f10Colours, add = TRUE, lwd = .5, pch = 19, cex = .7)
plot(AKpoly[1:100,], add = TRUE, border = 'black')
text(650000,1950000,'Steller Sea Lions', cex = 2.2, col = 'white')
text(650000,1710000,'Posterior', cex = 2.2, col = 'white')
text(650000,1450000,' ', cex = 2.2, col = 'white')
addBreakColorLegend(2000000, 1300000, 2200000, 2600000, f10$brks, 
    colors = pal, cex = 1.1, printFormat = "6.0")
    
xy = coordinates(trueCentroids)
#find 50 nearest neighbors for each hexagon
knnout = knn(xy, k = 50)
nb_50 = knnout$nn.idx
#create matrix to hold smoothed results
smoo50 = matrix(NA, nrow = dim(MCMCout$store_zmat)[1], 
	ncol = dim(MCMCout$store_zmat)[2])
#smooth over 50 nearest neighbors
start.time = Sys.time()
for(i in 1:dim(MCMCout$store_zmat)[1]) {
  for(j in 1:dim(MCMCout$store_zmat)[2]) {
    smoo50[i,j] = mean(temp[nb_50[i,],j])
  }
}
end.time = Sys.time()
difftime(end.time, start.time, units = 'mins')


# create objects from list items
MCMC_EJ_zmat = MCMCout$store_zmat
MCMC_EJ_sigma = MCMCout$store_sigma
MCMC_EJ_sigCAR = MCMCout$store_sigCAR
MCMC_EJ_rho =	MCMCout$store_rho
MCMC_EJ_beta0 =	MCMCout$store_beta0
MCMC_EJ_accept_samp = MCMCout$accept_samp
MCMC_EJ_accept_miss = MCMCout$accept_miss
MCMC_EJ_accept_beta0 = MCMCout$accept_beta0
MCMC_EJ_accept_sigma = MCMCout$accept_sigma
MCMC_EJ_accept_sigCAR = MCMCout$accept_sigCAR
MCMC_EJ_accept_rho = MCMCout$accept_rho
MCMC_EJ_smoo =	smoo50

#store the results in the package so they can be quickly accessed for graphics and 
#further processing
setwd(paste0('/home/jay/data/2019_papers/POP',
  '/POP_package/POP/data'))
save(MCMC_EJ_zmat, file = 'MCMC_EJ_zmat.rda')
save(MCMC_EJ_sigma, file = 'MCMC_EJ_sigma.rda')
save(MCMC_EJ_sigCAR, file = 'MCMC_EJ_sigCAR.rda')
save(MCMC_EJ_rho, file ='MCMC_EJ_rho.rda')
save(MCMC_EJ_beta0, file ='MCMC_EJ_beta0.rda')
save(MCMC_EJ_accept_samp, file = 'MCMC_EJ_accept_samp.rda')
save(MCMC_EJ_accept_miss, file = 'MCMC_EJ_accept_miss.rda')
save(MCMC_EJ_accept_beta0, file = 'MCMC_EJ_accept_beta0.rda')
save(MCMC_EJ_accept_sigma, file = 'MCMC_EJ_accept_sigma.rda')
save(MCMC_EJ_accept_sigCAR, file = 'MCMC_EJ_accept_sigCAR.rda')
save(MCMC_EJ_accept_rho, file = 'MCMC_EJ_accept_rho.rda')
save(MCMC_EJ_smoo, file = 'MCMC_EJ_smoo.rda')



