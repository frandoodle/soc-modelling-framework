combine_soc_depth_increments <- function(site_data_year) {
	# This function takes a group of site_data for one year with different ranges
	# of depth measurements and returns one row with the largest
	# possible depth range (e.g. 0-15cm + 15-30cm = 0-30cm)
	
	# Get SOC for each depth increment
	soc_by_depth <- site_data_year %>%
		select(soil_depth_min_cm, soil_depth_max_cm, soc_tha) %>%
		mutate(soc_per_half_cm = soc_tha / (soil_depth_max_cm-soil_depth_min_cm) / 2)
	# Early return is SOC is NA
	if(all(is.na(soc_by_depth$soc_tha))) {
		result <- site_data_year[1,] %>%
			mutate(soil_depth_min_cm = min(soc_by_depth$soil_depth_min_cm),
						 soil_depth_max_cm = max(soc_by_depth$soil_depth_max_cm),
						 soc_tha_lowest_depth = NA) %>%
			select(-soc_tha)
		return(result)
	}
	
	# Get the max measured depth
	max_depth <- max(soc_by_depth$soil_depth_max_cm)
	# Make a skeleton of SOC measurement depths by 0.5cm increments
	skeleton <- tibble(depth_min=seq(0,max_depth-0.5,0.5),
										 depth_max=seq(0.5,max_depth,0.5))
	# Merge SOC data into skeleton
	skeleton_filled <- skeleton
	for(row in 1:nrow(soc_by_depth)) {
		soc_by_depth_row <- soc_by_depth[row,]
		skeleton_filled[,row+2] <- ifelse(skeleton$depth_min >= soc_by_depth_row$soil_depth_min_cm & skeleton$depth_max <= soc_by_depth_row$soil_depth_max_cm,
																			soc_by_depth_row$soc_per_half_cm,
																			NA)
	}
	# Take the mean SOC from each depth increment
	soc_by_depth_increment_average <- skeleton_filled %>%
		mutate(average_soc_per_half_cm = rowMeans(skeleton_filled[3:length(skeleton_filled)], na.rm=T)) %>%
		select(depth_min,
					 depth_max,
					 average_soc_per_half_cm)
	
	soc_to_lowest_depth <- sum(soc_by_depth_increment_average$average_soc_per_half_cm)
	
	# Take the first row of site_data and replace every changed row to get result
	result <- site_data_year[1,] %>%
		mutate(lowest_depth_measured = max(soc_by_depth$soil_depth_max_cm),
					 soc_tha_lowest_depth = soc_to_lowest_depth) %>%
		select(-soc_tha, -bd_gcm3, -soil_depth_min_cm, -soil_depth_max_cm)
	
	return(result)
}
