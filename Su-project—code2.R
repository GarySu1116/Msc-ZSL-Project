# Load necessary libraries
install.packages("dplyr")
install.packages("ggplot")
install.packages("ggspatial")
install.packages("ggmap")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("patchwork")
library(patchwork)
library(readr)
library(dplyr)
library(vegan)
library(ggplot2)
library(sf)
library(ggspatial)
library(rlang)
library(terra)
library(ggmap)
library(rnaturalearth)
library(rnaturalearthdata)

# Read the dataset
# data <- read_csv("~/Desktop/UCL/Bios0034/coding/myproject/Code2")
data<-readxl::read_excel("Su-Data-MoreEnvironmentalVariables.xlsx")

# Extract species data (columns 9 to 23, excluding the "total" and "Density（total）" columns)
## I removed the column with only zeros
species_data_raw <- data[, 10:20]

## normalise species counts to a standard unit area (median distance is c.500m = 0.5 km so lets use that)
species_data <- 0.5 * species_data_raw / data$"Distance（km）"


# Calculate number of species per station (count non-zero entries)
data$species_count <- rowSums(species_data > 0)


# Calculate total abundance
data$total_abundance <- rowSums(species_data)

# Calculate Simpson diversity index
data$simpson_diversity <- diversity(species_data, index = "simpson")

# Calculate Shannon diversity index
data$shannon_diversity <- diversity(species_data, index = "shannon")

# Display the results
print(data)

#what are the most abundant VME taxa?
# Sum the abundance of each taxa across all stations
total_abundance_by_taxa <- colSums(species_data)

# Find the most abundant VME taxa
most_abundant_taxa <- sort(total_abundance_by_taxa, decreasing = TRUE)

# Display the results
print(most_abundant_taxa)

# do any taxa meet the benchmark abundance threshold of 1 per m2 at any station?
# Calculate the area covered at each station using the width of 1.49 meters
data$area_covered <- data$`Distance（km）` * 1000 * 1.49  # Convert km to meters and multiply by the width

# Calculate the abundance per m² for each taxa at each station

## now that the species_data set is normalised to 500m transect we need to do this differently

# abundance_per_m2 <- species_data / data$area_covered
abundance_per_m2 <-  species_data / 500

# Check if any taxa meet the benchmark abundance threshold of 1 per m² at any station
taxa_meeting_threshold <- abundance_per_m2 >= 1

# Find the taxa and stations where the threshold is met
taxa_and_stations_meeting_threshold <- which(taxa_meeting_threshold, arr.ind = TRUE)

# Display the results
if (nrow(taxa_and_stations_meeting_threshold) > 0) {
  result <- data.frame(
    Station = data$station[taxa_and_stations_meeting_threshold[, 1]],
    Taxa = colnames(species_data)[taxa_and_stations_meeting_threshold[, 2]],
    Abundance_per_m2 = abundance_per_m2[taxa_and_stations_meeting_threshold]
  )
  print(result)
} else {
  print("No taxa meet the benchmark abundance threshold of 1 per m² at any station.")
}

#do any stations have total VME density of above 1 per m2 (add all vme densities together)?
# Calculate total VME density for each station
data$total_vme_density <- rowSums(abundance_per_m2)

# Check if any stations have total VME density above 1 per m²
stations_meeting_threshold <- data$total_vme_density > 1

# Find the stations where the threshold is met
stations_above_threshold <- data[stations_meeting_threshold, c("station", "total_vme_density")]

# Display the results
if (nrow(stations_above_threshold) > 0) {
  print(stations_above_threshold)
} else {
  print("No stations have total VME density above 1 per m².")
}

#which stations have the most VME taxa present (count number of VME taxa present at each station)
# Calculate number of VME taxa present at each station (count non-zero entries)
data$species_count <- rowSums(species_data > 0)

# Identify the station(s) with the most VME taxa present
max_species_count <- max(data$species_count)
stations_with_most_vme_taxa <- data %>% filter(species_count == max_species_count)

# Display the results
print(paste("The maximum number of VME taxa present at any station is:", max_species_count))
print("Stations with the most VME taxa present:")
print(stations_with_most_vme_taxa[, c("station", "species_count")])


## when dealing with abundance we normally log transform the values to normalise the distribution
# Scatter plot of total abundance vs depth
ggplot(data, aes(x = Depth, y = log(total_abundance))) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Total Abundance vs Depth",
       x = "Depth (m)",
       y = "log(Total Abundance)") +
  theme_minimal()




## as with the scatter plot we should log transform our abundances
# Calculate total abundance
data$total_abundance <- rowSums(species_data)
data$log_abundance <- log(data$total_abundance + 1)

# Build the generalized linear model
## because your response is no longer straight counts the poisson family isn't right 
## also I've added fishing effort to the dataframe
# glm_total_abundance <- glm(log_abundance ~ Depth +Temperature + Salinity + Fishing , data = data, family = poisson())
glm_total_abundance <- glm(log_abundance ~ Depth + Temperature + Salinity + Fishing, data = data, family=gaussian())

# Summary of the model
summary(glm_total_abundance)


# Identify taxa that appear in 50% or more of the stations
taxa_presence <- colSums(species_data > 0) / nrow(species_data)
taxa_to_analyze <- names(taxa_presence[taxa_presence >= 0.5])


## change these as glm above

## new code - LOG TRANSFORM single species abundances
glm_results <- lapply(taxa_to_analyze, function(taxa) {
  formula <- as.formula(paste0("log(`", taxa, "`+1) ~ Depth + Temperature + Salinity + Fishing"))
  glm(formula, data = data, family = gaussian())
})


# Summarize the results
glm_summaries <- lapply(glm_results, summary)

# Display the summaries
for (i in 1:length(taxa_to_analyze)) {
  cat("\nSummary for", taxa_to_analyze[i], ":\n")
  print(glm_summaries[[i]])
}

###Figure 1 - map of area with simple background of depth contours and fishing areas, plotting the sites and displaying some information about diversity at these sites

# SIMPER analysis
# new code - group into fishing and not fishing
simper_result <- simper(species_data, data$Fishing>0)
summary(simper_result)

# Convert data to sf object
data_sf <- st_as_sf(data, coords = c("Start Long", "Start Lat"), crs = 4326)

# Transform to an appropriate projection system 
data_sf <- st_transform(data_sf, crs = "+proj=aeqd +lat_0=65 +lon_0=-57")

# Generate map
ggplot() +
  # Add site information
  geom_sf(data = data_sf, aes(size = total_abundance, color = simpson_diversity), shape = 21, fill = "red") +
  scale_size(range = c(2, 10)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.75, "in"),
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Map of Study Area with Depth Contours and Fishing Areas",
       x = "Longitude", y = "Latitude",
       size = "Total Abundance",
       color = "Simpson Diversity") +
  theme_minimal()

# Generate table with diversity and abundance information
summary_table <- data %>%
  select(station, Depth, total_abundance, simpson_diversity, shannon_diversity) %>%
  arrange(desc(total_abundance))

# Save table as CSV file
write.csv(summary_table, "diversity_summary_table.csv", row.names = FALSE)


# Load spatial data
bathy <- rast("Bathymetry.tif")
dcont <- vect("Contour100m.shp")
closed <- vect("GHL_trawling_closure_64.5N_to_68N_clipped below 600m_v2.shp")

# Transform data_sf to match the CRS of the other spatial data if needed
data_sf <- st_transform(data_sf, crs = st_crs(closed))

# Set up colour ramp
blue.col <- colorRampPalette(c("darkblue", "lightblue"))

# Plot background map limited to halibut area
plot(bathy, ylim = c(62, 66), xlim = c(-58, -54), col = blue.col(30), breaks = seq(-3000, 0, 100), legend = FALSE)

# Add contour lines
plot(dcont, add = TRUE, col = rgb(0, 0, 0, 0.5), lwd = 0.5)

# Add closed area
plot(closed, col = "darkgreen", add = TRUE, density = c(10, 20, 1))

# Now add your points
plot(data_sf, add = TRUE, pch = 21, cex = data$total_abundance / 450, col = "red", bg = rgb(1, 0, 0, alpha = 0.5))

# Add a legend for point sizes
legend("topright", legend = c(1000, 2000, 3000), pch = 21, pt.cex = c(10, 20, 30) / 100, col = "red", pt.bg = rgb(1, 0, 0, alpha = 0.5), title = "Total Abundance")

###Figure 2 MDS Analyst
#MDS
# Calculate the abundance per m² for each taxa at each station
abundance_per_m2 <- species_data / data$area_covered

# Calculate the average swept area
average_swept_area <- mean(data$area_covered)

# Standardize the taxa counts by multiplying by the average swept area
standardized_species_data <- abundance_per_m2 * average_swept_area

# Add the station information to the standardized table
standardized_table <- cbind(data$station, standardized_species_data)

# Export the standardized table to a CSV file
write.csv(standardized_table, "standardized_species_data.csv", row.names = FALSE)

# Read the standardized species data
standardized_data <- read_csv("standardized_species_data.csv")

# Remove the station column for MDS analysis
mds_data <- standardized_data[, -1]

# Perform MDS analysis
mds_result <- metaMDS(mds_data, distance = "bray", trymax = 100)

# Plot the MDS result
plot(mds_result, type = "t")

# Add environmental variables (e.g., Depth, Start Lat, Start Long) to the plot
env_data <- data[, c("Depth", "Temperature", "Salinity", "Fishing")]
envfit_result <- envfit(mds_result, env_data, permutations = 999)

# Plot the environmental vectors
plot(envfit_result, p.max = 0.05)

# Customize the MDS plot
plot(mds_result, type = "n")
points(mds_result, display = "sites", pch = 19, col = "blue")
plot(envfit_result, p.max = 0.05, col = "red")

# Add labels to the plot
ordilabel(mds_result, display = "sites", cex = 0.7, col = "black")




###Figure 3 - additional maps showing abundance of key taxa (e.g. Halipteris, Acanella, Flabellum) in the area.

# List of key taxa
key_taxa <- c("Asconema foliatum", "Acanella arbuscula", "Flabellum alabastrum","Anthoptilum grandiflorum")

# Generate maps for each key taxon
for (taxon in key_taxa) {
  data_sf[[taxon]] <- species_data[[taxon]]
  
  p <- ggplot() +
    geom_sf(data = data_sf, aes(size = !!sym(taxon), color = !!sym(taxon)), shape = 21, fill = "yellow") +
    scale_size(range = c(2, 10)) +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(1, "in"), pad_y = unit(1, "in"),
                           style = north_arrow_fancy_orienteering) +
    labs(title = paste("Map of", taxon, "Abundance in the Study Area"),
         x = "Longitude", y = "Latitude",
         size = "Abundance",
         color = "Abundance") +
    theme_minimal()
  
  ggsave(paste("Map_of_", gsub(" ", "_", taxon), "_Abundance.png", sep = ""), plot = p)
}


###2
# List of key taxa
key_taxa <- c("Asconema foliatum", "Acanella arbuscula", "Flabellum alabastrum", "Anthoptilum grandiflorum")

# Generate maps for each key taxon
for (taxon in key_taxa) {
  data_sf[[taxon]] <- species_data[[taxon]]
  
  p <- ggplot() +
    geom_sf(data = data_sf, aes(size = !!sym(taxon), color = !!sym(taxon)), shape = 21, fill = "red") +
    scale_size(range = c(2, 10)) +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(1, "in"), pad_y = unit(1, "in"),
                           style = north_arrow_fancy_orienteering) +
    labs(title = paste("Map of", taxon, "Abundance in the Study Area"),
         x = "Longitude", y = "Latitude",
         size = "Abundance",
         color = "Abundance") +
    theme_minimal() +
    
    # Use a simple longitude/latitude projection (WGS 84)
    coord_sf(crs = "+proj=longlat +datum=WGS84", expand = FALSE)
  
  ggsave(paste("Map_of_", gsub(" ", "_", taxon), "_Abundance.png", sep = ""), plot = p)
}





###Figure 4 - response plots for anything that shows up significant in your models (Based on previous result, I will consider Depth, Start Lat, and Start Long)


# Generate plot: Response of significant variables to total abundance
significant_variables <- c("Depth", "Temperature", "Salinity", "Fishing")  # List of significant variables

# Loop through each significant variable and generate plots
for (predictor in significant_variables) {
  if (predictor == "Fishing") {
    p <- ggplot(data, aes_string(x = paste0("log(", predictor, " + 1)"), y = "log(total_abundance)")) +
      geom_point() +
      geom_smooth(method = "glm", method.args = list(family = "gaussian"), col = "blue") +
      labs(title = paste("Response of Log Total Abundance to Log", predictor),
           x = paste("Log", predictor),
           y = "Log Total Abundance") +
      theme_minimal()
  } else {
    p <- ggplot(data, aes_string(x = paste0("`", predictor, "`"), y = "log(total_abundance)")) +
      geom_point() +
      geom_smooth(method = "glm", method.args = list(family = "gaussian"), col = "blue") +
      labs(title = paste("Response of Log Total Abundance to", predictor),
           x = predictor,
           y = "Log Total Abundance") +
      theme_minimal()
  }
  
  ggsave(paste("Response_Plot_Log_Total_Abundance_", gsub(" ", "_", predictor), ".png", sep = ""), plot = p)
}

###Table 1 - station metadata and taxon summary stats

# Extract station metadata (only keeping Depth, Start Lat, Start Long, and Fishing)
station_metadata <- data %>%
  select(station, Depth, `Temperature`, `Salinity`, Fishing)

# Calculate taxon summary stats
taxon_summary_stats <- species_data %>%
  summarise_all(list(
    mean = mean,
    sd = sd,
    min = min,
    max = max,
    median = median
  ))

# Combine station metadata and taxon summary stats
table_1 <- bind_cols(station_metadata, taxon_summary_stats)

# Display the table
print(table_1)

# Save as CSV file
write.csv(table_1, "Station_Metadata_and_Taxon_Summary_Stats.csv", row.names = FALSE)





station_metadata <- data %>%
  select(station, Depth, Temperature, Salinity, Fishing)


species_data_with_station <- bind_cols(station_metadata, species_data)


taxon_summary_stats <- species_data_with_station %>%
  group_by(station) %>%  
  summarise(across(everything(), list(
    mean = ~mean(.x, na.rm = TRUE),  
    sd = ~sd(.x, na.rm = TRUE),      
    min = ~min(.x, na.rm = TRUE),    
    max = ~max(.x, na.rm = TRUE),    
    median = ~median(.x, na.rm = TRUE) 
  ), .names = "{col}_{fn}"))  


table_1 <- left_join(station_metadata, taxon_summary_stats, by = "station")


print(table_1)

write.csv(table_1, "Station_Metadata_and_Taxon_Summary_Stats.csv", row.names = FALSE)


###Table 2 - summary of your GLMs showing estimate, error & p-values for each variable (1 row per model)

# List of species to analyze
species_list <- colnames(species_data)

# Create an empty dataframe to store results
glm_summary <- data.frame(Species = character(),
                          Variable = character(),
                          Estimate = numeric(),
                          Std_Error = numeric(),
                          P_Value = numeric(),
                          stringsAsFactors = FALSE)

# Fit GLM for each species and extract results
for (species in species_list) {
  formula <- as.formula(paste0("`", species, "` ~ Depth + Temperature + Salinity + Fishing"))
  model <- glm(formula, data = data, family = gaussian())
  summary_model <- summary(model)
  
  for (var in rownames(summary_model$coefficients)) {
    glm_summary <- rbind(glm_summary, data.frame(
      Species = species,
      Variable = var,
      Estimate = summary_model$coefficients[var, "Estimate"],
      Std_Error = summary_model$coefficients[var, "Std. Error"],
      P_Value = summary_model$coefficients[var, "Pr(>|t|)"]
    ))
  }
}

# View results
print(glm_summary)

# Save results to CSV file
write.csv(glm_summary, "GLM_Summary_Stats.csv", row.names = FALSE)
