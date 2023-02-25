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
