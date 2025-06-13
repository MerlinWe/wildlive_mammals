con2 <- con2 %>%
	filter(complete.cases(Species)) %>%
	mutate(Species = paste(Genus, Species, sep = " ")) %>%
	select(-Trivial)

unanimous <- con2 %>% 
	dplyr::select(External_ID, Species) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Species, 
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
batch1 <- subset(con2, External_ID %in% batch1) 
con2 <- subset(con2, !(External_ID %in% batch1$External_ID)) 

batch1 <- batch1 %>% 
	dplyr::select(External_ID, Link, Category, Family, Genus, Species,
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch1")

remainings <- con2 %>% 
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	dplyr::select(External_ID, Index, Species) %>% 
	pivot_wider(names_from = Index, 
							values_from = Species, 
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
colnames(batch2) <- c("External_ID", "Species")

id <- paste0(batch2$External_ID, batch2$Species) 
con2$ID <- paste0(con2$External_ID, con2$Species) 
batch2 <- distinct(subset(con2, ID %in% id), ID, .keep_all = TRUE)[ ,-15] 

batch2 <- batch2 %>%
	dplyr::select(External_ID, Link, Category, Family, Genus, Species,
								Station, Camera, DateTimeOriginal) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch2")

rm(id)
con2 <- bind_rows(batch1, batch2) 
no_consensus <- subset(remainings, percentage <= 70)

