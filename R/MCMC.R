#-------------------------------------------------------------------------------
#
#           MCMC
#
#-------------------------------------------------------------------------------

#' MCMC sampling for count regression with a constant mean, latent autocorrelated errors, and missing data
#'
#' MCMC sampling for count regression with a constant mean, latent autocorrelated errors, and missing data
#'
#' @param response vector with response variable.
#' @param datadist probability density function for data.  It must be parameterized with a mean parameter, and a variance parameter.  If there is no variance parameter, the argument must still be there, and just evaluated to NULL.
#' @param antilink function to change mean (linear model) on link scale to mean for datadist distribution. 
#' @param beta0ini initialize beta0
#' @param sigmaini initialize sigma for datadist, if it has one.
#' @param sigCARini initialize sigCAR
#' @param rhoini initialize rho. It must be one of the values in rhogrid
#' @param zini initialize latent autocorrelated errors z
#' @param rhogrid vector of rho values for lookup table
#' @param detsig vector of determinant values corresponding to rhogrid
#' @param SigiList list of inverse covariances, in sparse matrix form,
#         corresponding to rhogrid
#' @param eff_mat_offset matrix of samples from posterior distribution for effort.  Default is NULL, in which case no offset is used
#' @param ind_samp vector of indexes of sampled grid cells
#' @param ind_miss vector of indexes of unsampled grid cells
#' @param sample_sigma Logical.  Should sigma parameter be sampled (TRUE) or held at initial values (FALSE).
#' @param sample_sigCAR Logical.  Should sigCAR parameter be sampled (TRUE) or held at initial values (FALSE).
#' @param sample_rho Logical.  Should rho parameter be sampled (TRUE) or held at initial values (FALSE).
#' @param use_trunc Logical.  Should truncation for z (latent autocorrelated errors) be used.  If true, truncated to plus and minus the log of the largest observed value.
#' @param nMCMC number of MCMC samples to retain
#' @param n_z_iters number of latent autocorrelated errors (z) proposals per further MCMC updating
#' @param n_thin number of MCMC samples per retained MCMC sample
#' @param samp_sd tuning parameter for width of independent uniform proposal for Metropolis step of sampled z-values
#' @param miss_sd tuning parameter for width of independent uniform proposal for Metropolis step of missing z-values
#' @param beta0_sd tuning parameter for variance of independent normal proposal for Metropolis step of beta0
#' @param sigma_sd tuning parameter for variance of independent normal proposal for Metropolis step of sigma
#' @param sigCAR_sd tuning parameter for variance of independent normal proposal for Metropolis step of sigCAR
#' @param rho_indxpm tuning parameter for number of adjacent grid values, plus and minus, of proposal for Metropolis/Hastings step of rho
#'
#' @return a list of MCMC samples of the posteriors and acceptance rates of Metropolis proposals
#'
#' @author Jay Ver Hoef
#' @export 

MCMC = function(response, datadist = dpois, antilink = exp, 
	beta0ini, sigmaini = NULL, sigCARini, rhoini, zini, 
	rhogrid, detsig, SigiList, eff_mat_offset = NULL, ind_samp, ind_miss, 
	sample_sigma = FALSE, sample_sigCAR = TRUE,
	sample_rho = TRUE, use_trunc = FALSE, trunc_bound = NULL,
	nMCMC = 1000, n_z_iters = 10, n_thin = 1, samp_sd = .01, miss_sd = .01,
	beta0_sd = .1, sigma_sd = .1, sigCAR_sd = .02, rho_indxpm = 3)
{
	beta0 = beta0ini
	sigCAR = sigCARini
	sigma = sigmaini
	rho = rhoini
	z = zini
	Sigi = SigiList[[which(rhogrid == rho)]]
	
	if(use_trunc & is.null(trunc_bound)) {
		cat("\r", "Need non-NULL trunc_bound argument if use_trunc = TRUE")
		cat("\n")
		return('Need non-NULL trunc_bound argument if use_trunc = TRUE')
	}
	
	#-------------------------------------------------------------------------------
	# Do MCMC
	#-------------------------------------------------------------------------------

	#blank vectors/matrices/lists to store MCMC samples
	store_zmat = NULL
	store_sigCAR = NULL
	store_sigma = NULL
	store_rho = NULL
	store_beta0 = NULL
	accept_samp = NULL
	accept_miss = NULL
	accept_beta0 = NULL
	accept_sigma = NULL
	accept_sigCAR = NULL
	accept_rho = NULL
	
	start.time = Sys.time()
	for(kk in 1:nMCMC) {
		cat("\r", "MCMC iteration: ", kk)
		# thin the samples -- keep one out of every n_thin
		for(iter in 1:n_thin) {
      # if there is a matrix of MCMC samples of effort as an offset
			if(!is.null(eff_mat_offset)) {
					if(is.vector(eff_mat_offset)) rand_eff = eff_mat_offset[ind_samp]
					if(is.matrix(eff_mat_offset)) { 
						rand_eff = eff_mat_offset[ind_samp,
							sample(1:dim(eff_mat_offset)[2],1)]
					} 
			}
      #otherwise set effort as 1 (so it will be 0 on log scale)
			if(is.null(eff_mat_offset))	rand_eff = rep(1, times = length(ind_samp))
			# do n_z_iter samples of z for every other parameter sampling
			for(yiter in 1:n_z_iters) {
				# block MCMC for sampled y's using Metropolis
				# make a copy of z to modify as a proposal
				z.try = z
        # usual proposal from a normal distribution
        if(use_trunc == FALSE) {
          z.try[ind_samp] = rnorm(length(ind_samp), 
            mean = z[ind_samp], sd = samp_sd)
        }
				# uniform proposal if using truncated normal distribution for z
        if(use_trunc == TRUE) {
          lb = z[ind_samp] - samp_sd/2
          lb[z[ind_samp] - samp_sd/2 < -trunc_bound] = -trunc_bound
          lb[z[ind_samp] + samp_sd/2 > trunc_bound] = trunc_bound - samp_sd
          z.try[ind_samp] = lb + runif(length(ind_samp))*samp_sd 
        }
				# log Metropolis ratio
				# use uncertainty of effort surface
				LLdif = (sum(datadist(response[ind_samp], mu =
						antilink(beta0 + log(rand_eff) +
							z.try[ind_samp]), sigma = sigma, log = TRUE)) - 
							z.try %*% Sigi %*% z.try/(2*sigCAR^2)) - 
					(sum(datadist(response[ind_samp], mu =
						antilink(beta0 + log(rand_eff) +
							z[ind_samp]), sigma = sigma, log = TRUE)) - 
							z %*% Sigi %*% z/(2*sigCAR^2))
        #acceptance or not
        U <- log(runif(1))
				if(as.numeric(LLdif) > U) {
					z[ind_samp] <- z.try[ind_samp]
				}
				accept_samp = c(accept_samp, as.numeric(LLdif) > U)
				# block MCMC for unsampled y's using Metropolis	
				# make a copy of z to modify as a proposal
				z.try = z
       # usual proposal from a normal distribution
         if(use_trunc == FALSE) {
          z.try[ind_miss] = rnorm(length(ind_miss), 
            mean = z[ind_miss], sd = miss_sd)  
        }
				# uniform proposal if using truncated normal distribution for z
         if(use_trunc == TRUE) {
          lb = z[ind_miss] - miss_sd/2
          lb[z[ind_miss] - miss_sd/2 < -trunc_bound] = -trunc_bound
          lb[z[ind_miss] + miss_sd/2 > trunc_bound] = trunc_bound - miss_sd
          z.try[ind_miss] = lb + runif(length(ind_miss))*miss_sd 
        }
				LLdif = (-z.try %*% Sigi %*% z.try/(2*sigCAR^2)) + 
					(z %*% Sigi %*% z/(2*sigCAR^2))
       U <- log(runif(1))
				if(as.numeric(LLdif) > U) {
					z[ind_miss] <- z.try[ind_miss]
				}
			}
			accept_miss = c(accept_miss, as.numeric(LLdif) > U)
			# MCMC sample the intercept using Metropolis
			beta0.try = rnorm(1, beta0, beta0_sd)
			LLdif = sum(datadist(response[ind_samp], mu = antilink(beta0.try + 
					log(rand_eff) + z[ind_samp]), sigma = sigma,
					log = TRUE)) -
				sum(datadist(response[ind_samp], mu = antilink(beta0 + 
					log(rand_eff) + z[ind_samp]), sigma = sigma,
					log = TRUE))
			U <- log(runif(1))
			if(LLdif > U) {
				beta0 <- beta0.try
			}
			accept_beta0 = c(accept_beta0, as.numeric(LLdif) > U)
			# MCMC sample sigma using Metropolis
			if(sample_sigma) {
					sigma.try = exp(rnorm(1, log(sigma), sigma_sd))
				LLdif = sum(datadist(response[ind_samp], mu = antilink(beta0 + 
						log(rand_eff) + z[ind_samp]), sigma = sigma.try,
						log = TRUE)) -
					sum(datadist(response[ind_samp], mu = antilink(beta0 + 
						log(rand_eff) + z[ind_samp]), sigma = sigma,
						log = TRUE))
        U <- log(runif(1))
				if(LLdif > U) {
					sigma <- sigma.try
				}
				accept_sigma = c(accept_sigma, as.numeric(LLdif) > U)
			}
			if(sample_sigCAR) {
				# MCMC sample sigCAR using Metropolis
				logsigCAR.try = log(sigCAR) + rnorm(1, mean = 0, sd = sigCAR_sd)
				LLdif = (-z %*% Sigi %*% z/(2*exp(logsigCAR.try)^2) - length(z)*logsigCAR.try) + 
					(z %*% Sigi %*% z/(2*sigCAR^2) + length(z)*log(sigCAR))
        U <- log(runif(1))
        if(as.numeric(LLdif) > U) {
					sigCAR <- exp(logsigCAR.try)
				}
				accept_sigCAR = c(accept_sigCAR, as.numeric(LLdif) > U)
			}
			if(sample_rho) {
				# MCMC sample sigCAR using Hastings
				indx = which(rho == rhogrid)
				indx.pool = c(max(indx-rho_indxpm,1):max(indx-1,1),
					min(indx+1,length(rhogrid)):min(indx+rho_indxpm,length(rhogrid)))
				indx.try = sample(indx.pool,1)
				reverse.pool = c(max(indx.try-rho_indxpm,1):max(indx.try-1,1),
					min(indx.try+1,length(rhogrid)):min(indx.try+rho_indxpm,length(rhogrid)))
				rho.try = rhogrid[indx.try]
				LLdif = (-z %*% SigiList[[indx.try]] %*% z/(2*sigCAR^2) + 
					detsig[indx.try]/2) - 
					(-z %*% Sigi %*% z/(2*sigCAR^2) + 
					detsig[indx]/2)
				# Hastings uses the asymmetry in acceptance probabilities which
				# are 1/(number of indexes) [discrete uniform distribution]
				U <- log(runif(1))
				if(as.numeric(LLdif) + log(1/length(reverse.pool)) - 
						log(1/length(indx.pool)) > U) {
					rho = rho.try
					Sigi = SigiList[[indx.try]]
				}
				accept_rho = c(accept_rho, as.numeric(LLdif) + 
					log(1/length(reverse.pool)) - log(1/length(indx.pool)) > U)
			}	
		}
		#store results
		store_zmat = cbind(store_zmat,z)
		store_sigma = c(store_sigma, sigma)
		store_sigCAR = c(store_sigCAR, sigCAR)
		store_rho = c(store_rho, rho)
		store_beta0 = c(store_beta0, beta0)
	}
	cat("\n")
	list(store_zmat = store_zmat, store_sigCAR = store_sigCAR, 
		store_sigma = store_sigma, store_rho = store_rho, store_beta0 = store_beta0,
		accept_samp = accept_samp, 	accept_miss = accept_miss, 
		accept_beta0 = accept_beta0, accept_sigma = accept_sigma,
		accept_sigCAR = accept_sigCAR, accept_rho = accept_rho
	)
}
