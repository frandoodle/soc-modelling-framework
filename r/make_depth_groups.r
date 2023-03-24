make_depth_groups <- function(input) {
	
	z <- input %>%
		mutate(soil_depth = paste0(soil_depth_min_cm,"-",soil_depth_max_cm))
	
	number_of_measurements = z %>%
		filter(!is.na(soc_tha)) %>%
		nrow
	
	if(number_of_measurements == 0) {
		zero_combination_rows <- z %>%
			mutate(combination_group = 0,
						 lowest_depth = NA)
		return(zero_combination_rows)
	}
	
	
	one_combination_depths <- NULL
	two_combination_depths <- NULL
	three_combination_depths <- NULL
	
	one_combination_lowest_depth <- NULL
	two_combination_lowest_depth <- NULL
	three_combination_lowest_depth <- NULL
	
	if(nrow(z) >= 1) {
		one_combination <- combn(z$soil_depth, 1) %>%
			t %>%
			magrittr::set_colnames(c("V1")) %>%
			as_tibble(.name_repair = "unique") %>%
			separate(c("V1"), into = c("min","max"), sep = "-", convert = T, remove = F) %>%
			mutate(lowest_depth = (max) - (min)) %>%
			# Select the depth closest to 30
			ungroup() %>%
			mutate(lowest_depth_diff = abs(30-lowest_depth)) %>%
			filter(lowest_depth_diff == min(lowest_depth_diff))
		
		one_combination_lowest_depth <- one_combination$lowest_depth %>%
			unique()
		
		one_combination_depths <- one_combination %>%
			select(V1) %>%
			unlist
	}
	
	if(nrow(z) >= 2) {
		two_combination <- combn(z$soil_depth, 2) %>%
			t %>%
			magrittr::set_colnames(c("V1", "V2")) %>%
			as_tibble(.name_repair = "unique") %>%
			separate(c("V1"), into = c("min1","max1"), sep = "-", convert = T, remove = F) %>%
			separate(c("V2"), into = c("min2","max2"), sep = "-", convert = T, remove = F) %>%
			mutate(min_unique = min1 == min2,
						 max_unique = max1 == max2,
						 lowest_depth = (max1 + max2) - (min1 + min2)) %>%
			filter(!min_unique,
						 !max_unique) %>%
			# Select the depth closest to 30
			ungroup() %>%
			mutate(lowest_depth_diff = abs(30-lowest_depth)) %>%
			filter(lowest_depth_diff == min(lowest_depth_diff))
		
		two_combination_lowest_depth <- two_combination$lowest_depth %>%
			unique()
		
		two_combination_depths <- two_combination %>%
			select(V1, V2) %>%
			unlist
	}
	
	if(nrow(z) >= 3) {
		three_combination <- combn(z$soil_depth, 3) %>%
			t %>%
			magrittr::set_colnames(c("V1", "V2", "V3")) %>%
			as_tibble(.name_repair = "unique") %>%
			separate(c("V1"), into = c("min1","max1"), sep = "-", remove = F) %>%
			separate(c("V2"), into = c("min2","max2"), sep = "-", remove = F) %>%
			separate(c("V3"), into = c("min3","max3"), sep = "-", remove = F) %>%
			mutate(min_unique1 = min1 == min2,
						 min_unique2 = min1 == min3,
						 min_unique3 = min2 == min3,
						 max_unique1 = max1 == max2,
						 max_unique2 = max1 == max3,
						 max_unique3 = max2 == max3,
						 lowest_depth = (as.double(max1) + as.double(max2) + as.double(max3))
						 - (as.double(min1) + as.double(min2) + as.double(min3))) %>%
			filter(!min_unique1,
						 !min_unique2,
						 !min_unique3,
						 !max_unique1,
						 !max_unique2,
						 !max_unique3) %>%
			# Select the depth closest to 30
			ungroup() %>%
			mutate(lowest_depth_diff = abs(30-lowest_depth)) %>%
			filter(lowest_depth_diff == min(lowest_depth_diff))
		
		three_combination_lowest_depth <- three_combination$lowest_depth %>%
			unique()
		
		three_combination_depths <- three_combination %>%
			select(V1, V2, V3) %>%
			unlist
	}
	
	one_combination_rows <- z %>%
		filter(soil_depth %in% one_combination_depths) %>%
		mutate(combination_group = 1,
					 lowest_depth = one_combination_lowest_depth)
	two_combination_rows <- z %>%
		filter(soil_depth %in% two_combination_depths) %>%
		mutate(combination_group = 2,
					 lowest_depth = two_combination_lowest_depth)
	three_combination_rows <- z %>%
		filter(soil_depth %in% three_combination_depths) %>%
		mutate(combination_group = 3,
					 lowest_depth = three_combination_lowest_depth)
	
	result <- bind_rows(one_combination_rows, two_combination_rows, three_combination_rows) %>%
		filter(lowest_depth == max(lowest_depth))
	
	print(result$location_name)
	print(result$treatment_name)
	
	return(result)
}
