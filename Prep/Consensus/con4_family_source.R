con4 <- con4 %>%
	filter(complete.cases(Family)) %>%
	select(-Trivial, -Species, -Genus)

unanimous <- con4 %>% 
	dplyr::select(External_ID, Family) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Family, 
							names_prefix = "Class",
							values_fn = function(x) names(sort(table(x), decreasing = TRUE)[1]))

names <- colnames(unanimous) 
unanimous <- data.frame(t(apply(unanimous, 1, shift_to_left)))
colnames(unanimous) <- names 
unanimous <- unanimous %>% 
	filter(!rowSums(is.na(.[, 2:ncol(unanimous)])) == length(2:ncol(unanimous)))

unanimous_test <- unanimous %>%
	pmap_dfr(., ~ na.locf(c(...)) %>% as.list %>% as_tibble) %>%
	transform(Agree = apply(.[, 2:ncol(.)], 1,function(x) length(unique(x)) == 1)) %>%
	select(External_ID, Agree)
unanimous$Agree <- unanimous_test$Agree


batch1 <- subset(unanimous, Agree == TRUE)$External_ID 
batch1 <- subset(con4, External_ID %in% batch1) 
con4 <- subset(con4, !(External_ID %in% batch1$External_ID)) 

batch1 <- batch1 %>% 
	dplyr::select(External_ID, Link, Category, Family,
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch1")

remainings <- con4 %>% 
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	dplyr::select(External_ID, Index, Family) %>% 
	pivot_wider(names_from = Index, 
							values_from = Family, 
							names_prefix = "Class", 
							values_fn = function(x) names(sort(table(x), decreasing = TRUE)[1]))

External_ID <- remainings$External_ID
remainings <- apply(remainings, 1, function(row) {
	counts <- table(row)
	most_common <- names(which.max(counts))
	percentage <- counts[which.max(counts)] / sum(counts) * 100
	return(list(most_common = most_common, percentage = percentage))
})

remainings <- remainings %>%
	bind_rows() %>%  
	unnest_longer(percentage, indices_include = TRUE) %>% 
	mutate(External_ID = External_ID) %>%
	select(External_ID, most_common, percentage) %>%  
	mutate(most_common = as.character(most_common))

batch2 <- subset(remainings, percentage > 70)[, 1:2]
colnames(batch2) <- c("External_ID", "Family")
remainings <- subset(remainings, percentage < 70 )

id <- paste0(batch2$External_ID, batch2$Family) 
con4$ID <- paste0(con4$External_ID, con4$Family) 
batch2 <- distinct(subset(con4, ID %in% id), ID, .keep_all = TRUE)[ ,-15] 

batch2 <- batch2 %>%
	dplyr::select(External_ID, Link, Category, Family,
								Station, Camera, DateTimeOriginal) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch2")
rm(id)  

remainings <- con4 %>%
	filter(External_ID %in% remainings$External_ID) %>%
	filter(Review == 1) %>% 
	dplyr::select(External_ID, Family) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Family, 
							names_prefix = "Class",
							values_fn = function(x) names(sort(table(x), decreasing = TRUE)[1]))

review_test <- remainings %>%
	pmap_dfr(., ~ na.locf(c(...)) %>% as.list %>% as_tibble) %>%
	transform(Agree = apply(.[, 2:ncol(.)], 1,function(x) length(unique(x)) == 1)) %>%
	select(External_ID, Agree)
review <- remainings %>%
	left_join(review_test, by = "External_ID") %>%
	filter(Agree == TRUE)

batch3 <- filter(con4, External_ID %in% review$External_ID) %>%
	dplyr::select(External_ID, Link, Category, Family,
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch3")

con4 <- bind_rows(batch1, batch2, batch3) 

no_consensus <- remainings %>% 
	left_join(review_test, by = "External_ID") %>%
	filter(Agree == FALSE)