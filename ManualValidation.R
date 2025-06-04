##################################VALIDATION######################################
#################STEP 1: REMOVE DOUBLE COUNTS##############################
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(openxlsx)


#### Read in all BirdNET tables #####
folder_path <- "/Users/kelsey/Desktop/EppingBirdNETtables"
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {
  # Extract the file name without the extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Read the file into a data frame
  assign(file_name, read.delim(file_path))
}


#### Combine the dfs to determine total obs. and no. species #####
df_names <- ls(pattern = "BirdNET_SelectionTable")
total_table <- bind_rows(lapply(df_names, get))
nrow(total_table)
print(unique(total_table$Common.Name))
# Total obs. = 119720 and Total no. species = 236

##### If same species occurs within same WAV file, only keep one occurrence (with highest confidence) #####
folder_path <- "/Users/kelsey/Desktop/EppingBirdNETtables"
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {
  # Extract file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Remove 'BirdNETtable' part
  site_name <- sub("BirdNETtable", "", file_name)
  # Read file into data frame
  df <- read.delim(file_path)
  # Extract only the part of the 'Begin.Path' starting from '2025' onwards
  df$EndofBegin.Path <- sub(".*(2025.*)", "\\1", df$Begin.Path)
  # Filter so that if species is occurring multiple times within the same wav file (same site, day, and minute)
  # only the occurrence with the highest Confidence is kept
  df_filtered <- df %>%
    group_by(EndofBegin.Path, Species.Code) %>%
    slice_max(order_by = Confidence, n = 1) %>%
    ungroup()
  # Assign the filtered data frame to a variable with the file name
  assign(paste0(site_name, "_filtered"), df_filtered)
}


## Add a Site column to each df
filtered_df_names <- grep("_filtered$", ls(), value = TRUE)

for (df_name in filtered_df_names) {
  # Extract base name before '_filtered'
  base_name <- sub("_filtered$", "", df_name)
  # Retrieve df by name
  df <- get(df_name)
  # Add 'Site' column with the base name value
  df <- df %>%
    mutate(Site = base_name)
  # Assign df back to its original name
  assign(df_name, df)
}

# It made a dataframe called df but its just Weils so remove it
rm(df_filtered)





#### Combine filtered dfs #####
filtered_df_names <- grep("_filtered$", ls(), value = TRUE)
# Retrieve dfs
filtered_dfs <- mget(filtered_df_names)
# Combine
combined_df <- bind_rows(filtered_dfs)
# Export
write.table(combined_df, file = "/Users/kelsey/Desktop/epping_combined_filtered_data.txt", sep = "\t", row.names = FALSE)

# Total number of observations (each row is an observation)
nrow(combined_df)
# Number of unique species
length(unique(combined_df$`Common.Name`))

## Total obs. = 16423 and Total no. species = 236

##########STEP 2: SELECT FILES FOR VALIDATION VIA STRATIFIED SAMPLING################
setwd("/Users/kelsey/Desktop")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(stats)
library(lme4)
library(broom)
library(openxlsx)
library(sampling)


totalbirdnet_df <- read.delim("/Users/kelsey/Desktop/epping_combined_filtered_data.txt")


###### Create a table for all species occurring <10 times ######
under_ten <- totalbirdnet_df %>%
  group_by(Common.Name) %>%
  filter(n() < 10) %>%
  ungroup()

print(unique(under_ten$Common.Name))
nrow(under_ten)
#393 out of 16423 observations, 155 out of 236 species

under_ten_unique <- under_ten %>%
  distinct(Common.Name, .keep_all = TRUE)

under_ten_species <- under_ten_unique %>%
  dplyr::select(Common.Name)
#had to specify select from dplyr if not will clash with other packages

write.xlsx(under_ten_species, file = "/Users/kelsey/Desktop/under_ten_species.xlsx")






#### REMOVE SPECIES OCCURING <10 TIMES (PRE-VALIDATION STEP) ####
over_ten <- totalbirdnet_df %>%
  group_by(Common.Name) %>%
  filter(n() >= 10) %>%
  ungroup()

print(unique(over_ten$Common.Name))
nrow(over_ten)
#16030 observations and 81 species

####  obs across  species (after filtering for occur >10 times) ####
## Export this
write.csv(over_ten, file = "/Users/kelsey/Desktop/species_over_ten.csv")


confidence_summary <- over_ten %>%
  group_by(Common.Name) %>%
  summarise(mean_confidence = mean(Confidence, na.rm = TRUE),
            median_confidence = median(Confidence, na.rm = TRUE),
            sd_confidence = sd(Confidence, na.rm = TRUE),
            min_confidence = min(Confidence, na.rm = TRUE),
            max_confidence = max(Confidence, na.rm = TRUE),
            n = n()) %>%
  arrange(desc(mean_confidence))



#### SELECT OBS TO VALIDATE USING STRATIFIED RANDOM SAMPLING ####
# STRATA ARE BASED ON CONFIDENCE SCORE INTERVALS - VALIDATE MORE FROM HIGHER STRATA #

# Create function for stratified random sampling
stratified_sampling <- function(df, species, n_low = 2, n_mid = 3, n_high = 5) {
  df %>%
    arrange(Confidence) %>%
    mutate(
      Stratum = ntile(Confidence, 3),
      Stratum = case_when(
        Stratum == 1 ~ "Low",
        Stratum == 2 ~ "Mid",
        Stratum == 3 ~ "High"
      )
    ) %>%
    group_by(Stratum) %>%
    group_modify(~ sample_n(.x, size = min(n(), if_else(.y$Stratum == "Low", n_low, if_else(.y$Stratum == "Mid", n_mid, n_high)))))
}

# Apply to df
for_validation <- over_ten %>%
  group_by(Common.Name) %>%
  group_modify(~ stratified_sampling(.x, .y$Common.Name)) %>%
  ungroup()


#### Add a column for filename to allow extraction from HP Desktop folder of wav files ####
for_validation <- for_validation %>%
  mutate(filename = str_extract(Begin.Path, "[^\\\\]+$"))

table(for_validation$Common.Name)
#black redstart 8, grasshopper-warbler 9, collared-dove 9, eurasian curlew 9, gray heron 9, marsh tit 9, red kite 9, rock pigeon 8
table(for_validation$Stratum)
#393 high, 162 low, 243 mid

## Export file 
write.xlsx(for_validation, file = "/Users/kelsey/Desktop/files_to_validate.xlsx", rowNames = FALSE)




##### ADD COLUMN OF STRATUM MIN AND MAX CONFIDENCE VALUES #####

## editing stratified sampling to add columns in for_validation that have min and max
## confidence value intervals associated with each species' high mid and low strata
stratified_sampling <- function(df, species, n_low = 2, n_mid = 3, n_high = 5) {
  df <- df %>%
    arrange(Confidence) %>%
    mutate(
      Stratum = ntile(Confidence, 3),
      Stratum = case_when(
        Stratum == 1 ~ "Low",
        Stratum == 2 ~ "Mid",
        Stratum == 3 ~ "High"
      )
    )
  
  intervals <- df %>%
    group_by(Stratum) %>%
    summarize(
      Min_Confidence = min(Confidence),
      Max_Confidence = max(Confidence),
      .groups = 'drop'
    )
  
  df <- df %>%
    left_join(intervals, by = "Stratum")
  
  sampled_df <- df %>%
    group_by(Stratum) %>%
    group_modify(~ sample_n(.x, size = min(n(), if_else(.y$Stratum == "Low", n_low, if_else(.y$Stratum == "Mid", n_mid, n_high)))))
  
  return(sampled_df)
}

# Apply to df
for_validation <- over_ten %>%
  group_by(Common.Name) %>%
  group_modify(~ stratified_sampling(.x, .y$Common.Name)) %>%
  ungroup()

## Export file - overwriting one without min and max confidence interval columns
write.xlsx(for_validation, file = "/Users/kelsey/Desktop/files_to_validate.xlsx", rowNames = FALSE)





#### CREATE DF CONTAINING THE MIN CONFIDENCE VALUES ASSOCIATED WITH EACH SPECIES' HIGH MID LOW STRATA ####

# Define the stratified sampling function
stratified_sampling <- function(df) {
  df %>%
    arrange(Confidence) %>%
    mutate(
      Stratum = ntile(Confidence, 3),
      Stratum = case_when(
        Stratum == 1 ~ "Low",
        Stratum == 2 ~ "Mid",
        Stratum == 3 ~ "High"
      )
    )
}

# Apply the stratified sampling function to assign strata
stratified_df <- over_ten %>%
  group_by(Common.Name) %>%
  group_modify(~ stratified_sampling(.x)) %>%
  ungroup()

# Calculate the minimum confidence values for each stratum and species
strata_min_df <- stratified_df %>%
  group_by(Common.Name, Stratum) %>%
  summarize(Min_Confidence = min(Confidence), .groups = 'drop') %>%
  pivot_wider(names_from = Stratum, values_from = Min_Confidence, names_prefix = "min_") %>%
  rename(high_min = min_High, mid_min = min_Mid, low_min = min_Low)

# Export as csv
write.csv(strata_min_df, file = "/Users/kelsey/Desktop/min_conf_per_strata.csv")


#################STEP 3: PROCESS TABLE OF VALIDATED OBSERVATIONS######################################
setwd("/Users/kelsey/Desktop")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(stats)
library(lme4)
library(broom)
library(openxlsx)
library(sampling)

## Read in the validated files csv
validated <- read.csv("/Users/kelsey/Desktop/eppingvalidatedfiles.csv") #files marked positive, false or unsure
table(validated$TrueFalseUnsure)

## Conservative approach - convert all unsure into false
class(validated$TrueFalseUnsure)

conservative <- validated
conservative$TrueFalseUnsure[conservative$TrueFalseUnsure == "Unsure"] <- "FALSE"
table(conservative$TrueFalseUnsure)


## Find any species that were 100% inaccurate during validation
species_with_all_false <- conservative %>%
  group_by(Common.Name) %>%
  summarize(all_false = all(TrueFalseUnsure == FALSE)) %>%
  filter(all_false) %>%
  pull(Common.Name)
print(species_with_all_false)



#Read in dataset
birds_df <- read.csv("/Users/kelsey/Desktop/species_over_ten.csv")

## Remove ecologically improbable species
species_to_remove <- c("Akohekohe", "Black Woodpecker", "Brambling", "Common Crane", "Great Bittern", "Green-winged Teal", "Hooded Crow", "Northern Goshawk", "Short-toed Treecreeper", "Wood Warbler", "Common Redstart", "European Pied Flycatcher", "Redwing", "Wood Lark", "Melodious Warbler", "Black Redstart", "Eurasian Golden Oriole", "Common Redstart", "Eurasian Curlew", "Eurasian Woodcock", "Water Rail", "Hawfinch", "Little Ringed Plover", "Tree Pipit", "Marsh Tit", "European Turtle-Dove")

# Remove all occurrences of the specified species from the dataframe
birds_filtered <- birds_df %>%
  filter(!(Common.Name %in% species_to_remove))
print(unique(birds_filtered$Common.Name))
nrow(birds_filtered)
#### This leaves 15092 obs. across 56 species ####

#### Now perform precision calculations - and then filter remaining species ####

## Calculate precision as: p = Tp / (Tp + Fp)
## Using 10 BirdNET confidence score thresholds placed evenly between 0.5-0.95

# Remove ecologically improbable species from conservative df
species_to_remove <- c("Akohekohe", "Black Woodpecker", "Brambling", "Common Crane", "Great Bittern", "Green-winged Teal", "Hooded Crow", "Northern Goshawk", "Short-toed Treecreeper", "Wood Warbler", "Common Redstart", "European Pied Flycatcher", "Redwing", "Wood Lark", "Melodious Warbler", "Black Redstart", "Eurasian Golden Oriole", "Common Redstart", "Eurasian Curlew", "Eurasian Woodcock", "Water Rail", "Hawfinch", "Little Ringed Plover", "Tree Pipit", "Marsh Tit", "European Turtle-Dove")
conservative <- conservative %>%
  filter(!(Common.Name %in% species_to_remove))
print(unique(conservative$Common.Name))

# Define the confidence thresholds
confidence_thresholds <- seq(0.5, 0.95, by = 0.05)

# Function to calculate precision for a given species and confidence threshold
calculate_precision <- function(df, species, threshold) {
  filtered_df <- df %>%
    filter(Common.Name == species & Confidence >= threshold)
  
  if(nrow(filtered_df) == 0) {
    return(NA)
  }
  
  true_positives <- sum(filtered_df$TrueFalseUnsure == "TRUE")
  total_positives <- nrow(filtered_df)
  
  precision <- true_positives / total_positives
  return(precision)
}

# Get unique species names
species_list <- unique(conservative$Common.Name)

# Create empty list to store results
results <- list()

# Loop through each species and confidence threshold to calculate precision
for(species in species_list) {
  for(threshold in confidence_thresholds) {
    precision_value <- calculate_precision(conservative, species, threshold)
    results <- append(results, list(data.frame(
      Species = species,
      ConfidenceThreshold = threshold,
      Precision = precision_value
    )))
  }
}

# Combine the results into a df
precision <- do.call(rbind, results)
view(precision)
print(unique(precision$Species))

# Re-name Species column
colnames(precision)[colnames(precision) == "Species"] <- "Common.Name"



#### Find all species' minimum confidence threshold required for 90% precision ####
species_list <- unique(precision$Common.Name)

# Create empty dataframe to store results
min_confidence_thresholds <- data.frame(Common.Name = character(), MinThreshold = numeric(), stringsAsFactors = FALSE)

# Loop through each species and confidence threshold to calculate precision and find minimum threshold for 90% precision
for(species in species_list) {
  min_threshold <- NA
  for(threshold in confidence_thresholds) {
    precision_value <- calculate_precision(conservative, species, threshold)
    if (!is.na(precision_value) && precision_value >= 0.9) {
      min_threshold <- threshold
      break
    }
  }
  min_confidence_thresholds <- rbind(min_confidence_thresholds, data.frame(Common.Name = species, MinThreshold = min_threshold))
}

view(min_confidence_thresholds)

#Grasshopper warbler, Bullfinch, Gadwall, Garden warbler, Gray heron and rock pigeon are NA.
#none of the predictions reach 90% precision for these species

library(dplyr)

min_confidence_thresholds <- min_confidence_thresholds %>%
  filter(!is.na(MinThreshold))
view(min_confidence_thresholds)


####### FILTER BIRD DATA BASED ON MINIMUM CONFIDENCE THRESHOLDS PER SPECIES #######
# join the MinThreshold column to birds_filtered based on Common.Name
birds_filtered_with_thresholds <- birds_filtered %>%
  left_join(min_confidence_thresholds, by = "Common.Name")

# Filter rows in birds_filtered_with_thresholds where Confidence is above the MinThreshold
filtered_birds <- birds_filtered_with_thresholds %>%
  filter(Confidence >= MinThreshold)

view(filtered_birds)
table(filtered_birds$Common.Name)


#### EXPORT FINAL DATASET TO BE USED ####
write.csv(filtered_birds, file = "/Users/kelsey/Desktop/final_bird_dataset.csv")


# Create table of total species counts
species_count <- as.data.frame(table(filtered_birds$Common.Name))
species_count <- species_count[order(-species_count$Freq), ]
# Export
write.xlsx(species_count, "/Users/kelsey/Desktop/species_count_final.xlsx")
