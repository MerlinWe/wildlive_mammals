review <- con1 %>%
	filter(Review == 1) %>%
	dplyr::select(External_ID, Trivial) %>%
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Trivial, 
							names_prefix = "Class_")

review_test <- review %>%
	pmap_dfr(., ~ na.locf(c(...)) %>% as.list %>% as_tibble) %>%
	transform(Agree = apply(.[, 2:ncol(.)], 1,function(x) length(unique(x)) == 1)) %>%
	select(External_ID, Agree)

review <- review %>%
	left_join(review_test, by = "External_ID") %>%
	filter(Agree == TRUE)

batch1 <- subset(con1, External_ID %in% review$External_ID) 
con1 <- con1 %>%
	filter(!(External_ID %in% batch1$External_ID)) %>%
	mutate(Review = 0)

batch1 <-  batch1 %>% 
	dplyr::select(External_ID, Link, Category, Family, Genus, Species, Trivial, 
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch1")

no_review <- con1 %>% dplyr::select(External_ID, Trivial) %>% 
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	pivot_wider(names_from = Index, 
							values_from = Trivial, 
							names_prefix = "Class", 
							values_fn = function(x) names(sort(table(x), decreasing = TRUE)[1]))

shift_to_left <- function(row) {
	non_na_values <- row[!is.na(row)]
	shifted_row <- c(non_na_values, rep(NA, length(row) - length(non_na_values)))
	shifted_row }

names <- colnames(no_review) 
no_review <- data.frame(t(apply(no_review, 1, shift_to_left)))
colnames(no_review) <- names

no_review <- no_review %>% 
	filter(!rowSums(is.na(.[, 2:ncol(no_review)])) == length(2:ncol(no_review)))

no_review_test <-  no_review %>% 
	pmap_dfr(., ~ na.locf(c(...)) %>% 
					 	as.list %>% as_tibble)

no_review_test <- transform(no_review_test, Agree = apply(no_review_test[, 2:ncol(no_review_test)], 1, 
																													function(x) length(unique(x)) == 1)) 
no_review$Agree <- no_review_test$Agree

batch2 <- subset(no_review, Agree == TRUE)$External_ID 
batch2 <- subset(con1, External_ID %in% batch2) 
con1 <- subset(con1, !(External_ID %in% batch2$External_ID)) 

batch2 <- batch2 %>% 
	dplyr::select(External_ID, Link, Category, Family, Genus, Species, Trivial, 
								Station, Camera, DateTimeOriginal) %>%
	distinct(External_ID, .keep_all = TRUE) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch2")

remainings <- con1 %>% 
	group_by(External_ID) %>%
	mutate(Index = row_number()) %>%
	ungroup() %>%
	dplyr::select(External_ID, Index, Trivial) %>% 
	pivot_wider(names_from = Index, 
							values_from = Trivial, 
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

batch3 <- subset(remainings, percentage > 70)[, 1:2]
colnames(batch3) <- c("External_ID", "Trivial")

id <- paste0(batch3$External_ID, batch3$Trivial) 
con1$ID <- paste0(con1$External_ID, con1$Trivial)
batch3 <- distinct(subset(con1, ID %in% id), ID, .keep_all = TRUE)[ ,-15]

batch3 <- batch3 %>%
	dplyr::select(External_ID, Link, Category, Family, Genus, Species, Trivial, 
								Station, Camera, DateTimeOriginal) %>%
	mutate(DateTimeOriginal = as.POSIXct(DateTimeOriginal, format = "%Y:%m:%d %H:%M:%S")) %>%
	mutate(Batch = "batch3")

rm(id)

con1 <- bind_rows(batch1, batch2, batch3) 
no_consensus <- subset(remainings, percentage <= 70)

