con3 <- con3 %>%
	filter(complete.cases(Genus)) %>%
	select(-Trivial, -Species)

review_weighted <- con3 %>%
	filter(Review == 1) %>%
	group_by(External_ID) %>%
	summarise(expand = n()) %>%
	filter(expand < 11) %>%
	right_join(filter(con3, Review == 1), by = "External_ID") %>%
	mutate(expand = ceiling(expand*.3)) %>%
	filter(complete.cases(expand)) %>% 
	uncount(weights = expand) 

con3_weighted <- bind_rows(review_weighted, filter(con3, Review != 1))
rm(review_weighted) 

unanimous <- con3_weighted %>% 
	dplyr::select(External_ID, Genus) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Genus, 
							names_prefix = "Class",
							values_fn = function(x) names(sort(table(x), decreasing = TRUE)[1]))

names <- colnames(unanimous) 
unanimous <- data.frame(t(apply(unanimous, 1, shift_to_left))) 
colnames(unanimous) <- names 
unanimous <- unanimous %>% 
	filter(!rowSums(is.na(.[, 2:ncol(unanimous)])) == length(2:ncol(unanimous))) %>% 
	select(where(~!all(is.na(.))))

unanimous_test <- unanimous %>%
	pmap_dfr(., ~ na.locf(c(...)) %>% as.list %>% as_tibble) %>%
	transform(Agree = apply(.[, 2:ncol(.)], 1,function(x) length(unique(x)) == 1)) %>%
	select(External_ID, Agree)

unanimous$Agree <- unanimous_test$Agree
batch1 <- subset(unanimous, Agree == TRUE)$External_ID 
batch1 <- subset(con3, External_ID %in% batch1) 
con3_weighted <- subset(con3_weighted, !(External_ID %in% batch1$External_ID))

batch1 <- batch1 %>% 
	dplyr::select(External_ID, Link, Category, Family, Genus, 
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch1")

remainings <- con3_weighted %>% 
	filter(!(External_ID %in% batch1$External_ID)) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	dplyr::select(External_ID, Index, Genus) %>% 
	pivot_wider(names_from = Index, 
							values_from = Genus, 
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
colnames(batch2) <- c("External_ID", "Genus")

id <- paste0(batch2$External_ID, batch2$Genus) 
con3$ID <- paste0(con3$External_ID, con3$Genus) 
batch2 <- distinct(subset(con3, ID %in% id), ID, .keep_all = TRUE)[ ,-15] 

batch2 <- batch2 %>%
	dplyr::select(External_ID, Link, Category, Family, Genus,
								Station, Camera, DateTimeOriginal) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch2")

rm(id)

con3 <- bind_rows(batch1, batch2) 
no_consensus <- subset(remainings, percentage <= 70)
