validation <- function(site_data,
											 climate_data,
											 initial_c,
											 distribution,
											 model = "ipcct2")
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
	
	result_list <- list()
	
	for(site_n in 1:length(site_data)) {
		validation_results_list <- list()
		site_identifier <- site_data[[site_n]] %>%
			select(Exp_ID, location_name, treatment_name, treatment_number) %>%
			unique
		
		ncores=parallel::detectCores()-2
		cl=parallel::makeCluster(ncores)
		doParallel::registerDoParallel(cl)
		
		model_results=foreach(parameters=1:nrow(distribution),
													.packages = c("parallel", 
																				"doParallel", 
																				"tidyverse"),
													.export = c("run_ipcct2",
																			"IPCCTier2SOMmodel")) %dopar%
			
			do.call("run_ipcct2",
							append(list(site_data[[site_n]],
													climate_data,
													init_active = initial_c[[site_n]]$init_active,
													init_slow = initial_c[[site_n]]$init_slow,
													init_passive = initial_c[[site_n]]$init_passive,
													return_site_inputs = TRUE),
										 parameters))
		
		valid <- model_results %>%
			purrr::map(function(x)
			{
				a <- x %>%
					select(year, soc_total, soc_tha_30cm) %>%
					na.omit()
				b <- validation_calculate_stats(simulated = a$soc_total/10,
																		 measured = a$soc_tha_30cm)
				return(b)
			}) %>%
			bind_rows
		
		validation_results_list[[site_n]] <- site_identifier %>%
			bind_cols(valid)
		stopCluster(cl)
		
		
		validation_stats <- validation_results_list %>%
			bind_rows
		result_list[[site_n]] <- validation_stats
	}
	
	result_list %>%
		bind_rows %>%
		return()
}

validation_calculate_stats <- function(simulated, measured) {
	# Takes: Two vectors
	# Returns: Tibble of all validation stats including root mean square error (RMSE),
	# 				 normalized RMSE (NRMSE), and Nash-Sutcliffe Efficiency (NSE)
	
	# RMSE
	P <- simulated
	O <- measured
	PminusO <- P - O
	PminusOsquared <- PminusO^2
	PminusOsquaredsums <- sum(PminusOsquared)
	n <- length(PminusO)
	RMSE <- sqrt((PminusOsquaredsums)/n)
	
	# NRMSE
	Obar <- mean(O)
	NRMSE <- RMSE / Obar
	
	# NSE
	OiminusObar <- O - Obar
	OiminusObarsquared <- OiminusObar^2
	OiminusObarsquaredsums <- sum(OiminusObarsquared, na.rm = TRUE)
	NSE <- ((OiminusObarsquaredsums) - (PminusOsquaredsums))/OiminusObarsquaredsums
	
	all_stats <- c(RMSE = RMSE,
								 NRMSE = NRMSE,
								 NSE = NSE)
	
	return(all_stats)
}
