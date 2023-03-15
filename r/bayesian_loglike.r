loglik=function(m,o){
  if(length(m)!=length(o)){
     print("Inequal number of modeled and observed values, cannot proceed")
     return()
   }
  
  res=log(m)-log(o)
  sigma=sqrt(mean(res^2))
  n=length(m)
  lk=-n*log(sigma)-(1/(2*sigma^2))*sum(res^2)
  return(lk)
  
}

run_ipcct2_calculate_loglik <- function(site_data,
																				climate_data,
																				init_active,
																				init_slow,
																				init_passive,
																				parameters) 
{
	id <- parameters[[1]] #gets id, hopefully (I think foreach::foreach coerces rows into unnamed vectors)
	modelled <- do.call("run_ipcct2",
											append(list(site_data = site_data,
																	climate_data = climate_data,
																	init_active = init_active,
																	init_slow = init_slow,
																	init_passive = init_passive),
														 parameters))
	actuals <-  site_data %>%
		mutate(POLYID = as.character(POLYID)) %>%
		select(site = POLYID, year = year_name,  actual = soc_tha_30cm)
	
	model_actual <- modelled %>%
		full_join(actuals, by=c("site", "year")) %>%
		select(site, year, soc_total, actual) %>%
		filter(!is.na(actual))
	
	loglike <- loglik(model_actual$soc_total, model_actual$actual)
	
	output <- list(model_return = bind_cols(modelled, parameters, site_data),
								 loglik = tibble(id, loglik = loglike))
	
	return(output)
}

run_loglike_parallel <- function(site_data,
																 climate_data,
																 initial_c,
																 parameter_sample)
{
	# site_data can either be a data.frame, or a list of data.frames
	if(!inherits(site_data, "list")) {
		site_data <- list(site_data)
	}
	#initial_c should always be a list of lists
	if(!any(sapply(initial_c, is.list))) {
		initial_c <- list(initial_c)
	}
	if(length(site_data) != length(initial_c)) {
		stop("site_data and initial_c cannot be different lengths")
	}
	
	#=================================================================================
	# Run the model and calculate log-likelihood
	#    - likelihoods were calculated assuming that the error (modeled - mseasured) are iid 
	#=================================================================================
	likelihood <- NULL
	likelihood_list <- list()
	model_return_list <- list()
	for(site_n in 1:length(site_data)) {
		# Begin parallel
		ncores=parallel::detectCores()-2
		cluster=parallel::makeCluster(ncores)
		doParallel::registerDoParallel(cluster)
		
		likelihood=foreach(i=1:nrow(parameter_sample), 
											 #.combine = rbind, 
											 .packages = c("parallel", 
											 							"doParallel", 
											 							"tidyverse"),
											 .export = c("run_ipcct2",
											 						"IPCCTier2SOMmodel",
											 						"loglik",
											 						"run_ipcct2_calculate_loglik")) %dopar%
			
			run_ipcct2_calculate_loglik(site_data = site_data[[site_n]],
																	climate_data = climate_data,
																	init_active = initial_c[[site_n]]$init_active,
																	init_slow = initial_c[[site_n]]$init_slow,
																	init_passive = initial_c[[site_n]]$init_passive,
																	parameters = parameter_sample[i,])
		stopCluster(cluster)
		# End parallel
		likelihood_list[[site_n]] <- likelihood %>%
			purrr::map(~.$loglik) %>%
			bind_rows
		model_return_list[[site_n]] <- likelihood %>%
			purrr::map(~.$model_return) %>%
			bind_rows
		# Status
		print(paste0("site ", site_n, "/", length(site_data), " (sample_size = ",nrow(parameter_sample),")"))
	}
	
	return(list(likelihood = likelihood_list,
							model_return = model_return_list))
}