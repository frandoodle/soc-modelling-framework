validation <- function(site_data,
											 climate_data,
											 initial_c,
											 model)
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
	
	for(site_n in 1:length(site_data)) {
		result <- run_ipcct2(site_data = site_data[[site_n]],
												 climate_data = climate_data,
												 init_active = initial_c[[site_n]]$init_active,
												 init_slow = initial_c[[site_n]]$init_slow,
												 init_passive = initial_c[[site_n]]$init_passive,
												 return_site_inputs = TRUE) %>%
			select(year, soc_total, soc_tha_30cm) %>%
			na.omit()
	}
	
	validation_stats <- validation_calculate_stats(simulated = result$soc_total,
																								 measured = result$soc_tha_30cm)
	
	return(validation_stats)
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
