source(here::here("r/bayesian_validation_calculate_stats.r"))

validation <- function(model_return)
{
	
	y <- model_return %>%
		group_by(SampleID) %>%
		select(SampleID, Exp_ID, location_name, treatment_name, treatment_number, year, soc_total, soc_tha_30cm) %>%
		na.omit() %>%
		group_split
	
	yy <- y %>%
		purrr::map(function(data) {
			data %>%
				mutate(soc_total_delta_cumulative = soc_total - first(soc_total),
							 soc_tha_30cm_delta_cumulative = soc_tha_30cm - first(soc_tha_30cm),
							 
							 soc_total_delta_progressive = soc_total - lag(soc_total),
							 soc_tha_30cm_delta_progressive = soc_tha_30cm - lag(soc_tha_30cm)) %>%
				# Set the first value for each delta soc column to NA so they aren't included in validation
		mutate(soc_total_delta_cumulative = ifelse(row_number() == 1, NA,soc_total_delta_cumulative),
					 soc_tha_30cm_delta_cumulative = ifelse(row_number() == 1, NA,soc_tha_30cm_delta_cumulative),
					 soc_total_delta_progressive = ifelse(row_number() == 1, NA,soc_total_delta_progressive),
					 soc_tha_30cm_delta_progressive = ifelse(row_number() == 1, NA,soc_tha_30cm_delta_progressive))
		})
	yyy <- yy %>%
		purrr::map(function(data) {
			identifier <- data %>%
				select(SampleID, Exp_ID, location_name, treatment_name, treatment_number) %>%
				unique
			
			validation_stocks <- validation_calculate_stats(simulated = data[["soc_total"]]/10,
																											measured = data[["soc_tha_30cm"]]) %>%
				enframe %>%
				pivot_wider(names_from = name) %>%
				rename_with(~ paste0(., "_stocks"))
			
			validation_delta_cumulative <- validation_calculate_stats(simulated = data[["soc_total_delta_cumulative"]]/10,
																																measured = data[["soc_tha_30cm_delta_cumulative"]]) %>%
				enframe %>%
				pivot_wider(names_from = name) %>%
				rename_with(~ paste0(., "_delta_cumulative"))
			validation_delta_progressive <- validation_calculate_stats(simulated = data[["soc_total_delta_progressive"]]/10,
																																 measured = data[["soc_tha_30cm_delta_progressive"]]) %>%
				enframe %>%
				pivot_wider(names_from = name) %>%
				rename_with(~ paste0(., "_delta_progressive"))
			
			return(bind_cols(identifier, validation_stocks, validation_delta_cumulative, validation_delta_progressive))
		}) %>%
		bind_rows
	
	return(list(stocks = yy,
							validation = yyy))
}
