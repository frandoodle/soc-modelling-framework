#=================================================================================
# Sampling Importance Resampling algorithm with Latin hypercube sampling in R
#=================================================================================
library(lhs)
library(doParallel)
library(parallel)
library(foreach)
library(here)

library(ggplot2)

source(here::here("r/ipcct2_run.r"))
source(here::here("r/bayesian_loglike.r"))

sir <- function(site_data,
								climate_data,
								initial_c,
								parameter_bounds,
								sample_size = 100,
								resample_size = 10) {
	#=================================================================================
	# read prior distribution
	# (Required columns: Parameter, value, lower, upper)
	paramBounds <- parameter_bounds
	#=================================================================================
	# names of parameters that are allowed to vary
	varSI       <- paramBounds$Parameter
	nParams     <- length(varSI)
	#=================================================================================
	# LHS sampling for SIR
	#=================================================================================
	# sample size (1000000 samples were used in Gurung et al., 2020)
	n <- sample_size
	X1 <- randomLHS(n = n, k = nParams)
	# transform standard uniform LHS to prior distribution 
	Y1 <- matrix(0, nrow = nrow(X1), ncol = nParams)
	for(i in 1:nParams){
		pos <- which(paramBounds$Parameter == varSI[i])
		lower <- paramBounds[pos, "lower"]
		upper <- paramBounds[pos, "upper"]
		Y1[, i] <- qunif(X1[, i], min = lower, max = upper)
	}
	X <- as.data.frame(Y1)
	names(X) <- varSI
	X <- cbind("SampleID" = 1:nrow(X), X)
	
	loglike_return <- run_loglike_parallel(site_data = calibration_data,
											 climate_data = climate_data,
											 initial_c = initial_c_calibration_spinup,
											 parameter_sample = X)
	
	model_return_list <- loglike_return$model_return
	
	Lkhood1 <- loglike_return$likelihood %>%
		bind_rows %>%
		# Remove invalid values of log-likelihood
		mutate(loglik = ifelse(loglik == -Inf, NA, loglik)) %>%
		group_by(id) %>%
		summarise(loglik = mean(loglik, na.rm=T)) %>%
		mutate(weights = exp(loglik)/sum(exp(loglik)))
	
	#=================================================================================
	# sample without replacement (a resampling of 1000 were used in Gurung et al., 2020)
	nsamp <- resample_size 
	sampIndx <- sample(1:nrow(Lkhood1), 
										 size = nsamp,
										 replace = FALSE,
										 prob = Lkhood1$weights)
	PostTheta <- as.data.frame(X[sampIndx,])
	
	model_return_list_posterior <- model_return_list %>%
		purrr::map(~filter(., SampleID %in% PostTheta$SampleID))
	
	return(list(prior = X, posterior = PostTheta,
							model_return_prior = model_return_list,
							model_return_posterior = model_return_list_posterior))
}