# This script runs through every operation that is done to turn site_data_raw
# into site_data_processed

process_site_data <- function(site_data_raw) {
	
	site_data_processed_1 <- site_data_raw %>%
		select(-treatment_name) %>%
		# Default tillage is conventional till
		mutate(tillage = ifelse(is.na(tillage), "CT", tillage)) %>%
		# Remove sites where there are no measured values, or only 1 measured value
		group_by(TrtID, replication_number, block) %>%
		filter(!all(is.na(soc_tha))) %>%
		filter(!sum(!is.na(soc_tha), na.rm=T) <= 1)
	
	
	# Aggregate data about multiple crops
	# (If there are multiple fields where each sequence of the rotation is
	# being grown for the same year, we get the average yield/residue across
	# each field/crop. If crops are grown on the SAME field,
	# then we get the sum (e.g. cover crops))
	
	# These are all the columns that we use for summarising intercropping / rotation
	# related data. We will either be taking the mean value or the first value of these
	# columns per treatment per year.
	current_crop_columns = c("crop",
													 "ICASA_code",
													 "AAFC_code")
	
	yield_residue_columns <- c("grain_yield_kgha",
														 "crop_residue_kgha",
														 "hay_yield_kgha",
														 "hay_residue_kgha",
														 "root_biomass_kgha",
														 "straw_biomass_kgha",
														 "silage_yld_kgha",
														 #"biomass_returned_kgha",
														 "crop_residue_avg_kgha",
														 "above_ground_residue_kgha",
														 "roots_residue_kgha",
														 "root_exudates_kgha",
														 "manure_kgha")
	
	fertilizer_columns <- c("N_fertilizer_kgha",
													"N_fertilizer_kgha_2",
													"N_fertilizer_kgha_3",
													"P_fertilizer_kgha",
													"P_fertilizer_kgha_2",
													"P_fertilizer_kgha_3",
													"K_fertilizer_kgha",
													"K_fertilizer_kgha_2",
													"S_fertilizer_kgha",
													"S_fertilizer_kgha_2")
	
	site_data_processed_2 <- site_data_processed_1 %>%
		# 1. Summarise intercropping (same field)
		group_by_all %>%
		ungroup(!!current_crop_columns, !!yield_residue_columns, !!fertilizer_columns) %>%
		summarise(across(current_crop_columns, first),
							across(yield_residue_columns, sum),
							across(fertilizer_columns, sum)) %>%
		# 2. Summarise rotation sequences
		group_by_all %>%
		ungroup(plot_id, !!current_crop_columns, !!yield_residue_columns, !!fertilizer_columns) %>%
		summarise(across(current_crop_columns, first),
							across(yield_residue_columns, mean),
							across(fertilizer_columns, mean))
	
	
	# Remove all sites that have too many rows per treatment per year.
	# Select sites where there are over 1 row per trmt per year (per replication per block per sampling depth)
	offending_sites <- site_data_processed_2 %>%
		group_by(TrtID, replication_number, block, year_name, soil_depth_min_cm, soil_depth_max_cm, soc_tha) %>%
		tally() %>% #Gets the number of rows in the group
		filter(n > 1) %>%
		pull(TrtID) %>%
		unique()
	
	warning(paste0("Removing sites ",paste0(offending_sites, collapse=", "), " due too many rows per treatment per year"))
	
	site_data_processed_3 <- site_data_processed_2 %>%
		group_by(TrtID, replication_number, block, year_name) %>%
		filter(!TrtID %in% offending_sites)
	
	## Get measured SOC to lowest depth, then correct for 30cm
	
	# 1. Get the sum of measurement depths to get 0-lowest depth
	site_data_processed_4 <- site_data_processed_3 %>%
		group_by(TrtID, replication_number, block, year_name) %>%
		group_split %>%
		purrr::map(combine_depth_increments) %>%
		bind_rows
	
	# 2. Correct to 30cm
	site_data_processed_5 <- site_data_processed_4 %>%
		group_by(TrtID, replication_number, block) %>%
		mutate(lowest_depth_difference = 30 - lowest_depth_measured,
					 first_soc_treatment = first(na.omit(soc_tha_lowest_depth)),
					 soc_dsm_difference = initial_dsm - first_soc_treatment,
					 soc_tha_30cm = ifelse(lowest_depth_difference != 0,
					 											soc_tha_lowest_depth + soc_dsm_difference,
					 											soc_tha_lowest_depth))
	
}
