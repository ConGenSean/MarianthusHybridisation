####Install the necessary packages####
###You only need to do this once - once installed, they can be loaded using 'library'###
install.packages("terra", dependencies = T) #The main spatial analysis package
install.packages("sf", dependencies = T) #One of the core spatial packages
install.packages("galah", dependencies = T) #The Atlas of Living Australia package
install.packages("tidyverse", dependencies = T) #The tidyverse contains multiple packages designed for data wrangling in R
install.packages("reshape2", dependencies = T) #Another package for managing matrices and data in R

####Load relevant packages####
library(terra)
library(sf)
library(galah)
library(tidyverse)
library(reshape2)
#If at any stage in your code you get an error which says "could not find function", it means you do not have the necessary package loaded
#Google which package it is a part of, and then install and load the library as above if needed

####Set your working directory####
###This will define where R loads and exports files from###
setwd() #This sets the working directory to save outputs to - change to wherever you choose

####Script to import Atlas of Living Australia data####
###This code will take a list of species and collect the GPS data from ALA###
galah_config(email = "your.email.here@our.ecu.edu.au") #Connect to your ALA account

#Create a list of all the Marianthus species on ALA - we will 'loop' over each one to download the records
species_list <- c("Marianthus bignoniaceus", "Marianthus bicolor", "Marianthus erubescens","Marianthus candidus", "Marianthus drummondianus",
                  "Marianthus tenuis", "Marianthus coeruleopunctatus", "Marianthus sylvaticus","Marianthus ringens", "Marianthus granulatus",
                  "Marianthus microphyllus", "Marianthus mollis", "Marianthus paralius","Marianthus aquilonaris", "Marianthus dryandra")

####For each species in the list, downloading the ALA data####
for(species in species_list){
  cat("Downloading data for ",species,".\n", sep = "")
  init <- galah_call() %>%
    galah_identify(species) %>% #Select species
    atlas_occurrences() %>% #Download the occurrence data
    drop_na() #Remove those missing coordinates
  occ <- data.frame(init$decimalLatitude, 
                    init$decimalLongitude, 
                    init$scientificName) #Only keep the relevant columns
  
  cat("Initial dataset: ", nrow(occ)," observations.\n", sep = "")
  #These steps identify outliers and removes them
  #Specifically, it removes coordinates which are more than 3 standard deviations away from the rest of the data
  #May be tricky with different sampling so good to verify afterwards
  lat_scale <- scale(occ$init.decimalLatitude) #Rescale latitude to mean +/- standard deviations
  lat_outliers <- occ[abs(lat_scale) > 3,] #Check for outliers
  lon_scale <- scale(occ$init.decimalLongitude) #Rescale longitude to mean +/- standard deviations
  lon_outliers <- occ[abs(lon_scale) > 3,] #Check for outliers
  occ_no_out <- occ %>%
    filter(!init.decimalLatitude %in% lat_outliers$init.decimalLatitude) %>% #Remove latitude outliers
    filter(!init.decimalLongitude %in% lon_outliers$init.decimalLongitude) #Remove longitude outliers
  cat("Removing outliers: ", nrow(occ)-nrow(occ_no_out), " outliers found.\n", sep = "")
  
  colnames(occ_no_out) <- c("Latitude", "Longitude", "Species")
  #Save the occurrence data as a comma-delimited file for reference - can be opened in Excel
  spec_name <-gsub(pattern = " ", replacement = "_", x = species) #Replace the space in the species name with an underscore
  write.csv(x = occ_no_out, file = paste0(spec_name,"_ALA_occurrence_data.csv"), quote = F, row.names = F)
}

####Spatial analysis####
#First, list your input files - these should be saved in the same directory as your working directory above
#Make sure you have not changed their names or locations - if you have added FloraBase data to those files, that is okay
poly_list <- list.files(path = "./", pattern = "ALA_occurrence_data.csv", full.names = T)

buffer_width <- 2000 #Set the area of buffer around the point you would like to consider here
#Current value = 2km: this was generated based on tests with Marianthus aquilonaris
#If some of your files look weird, feel free to change this number

poly_df = NULL #Create an empty object to store the polygon results in

###For each species, create a polygon using a buffer radius###
for(species in poly_list){
  spec_split <- str_split(string = species, pattern = "_")[[1]] #Break out the elements of the occurrence data filename
  spec_name <- paste0(gsub(pattern = "./", replacement = "", x = spec_split[1]),"_",spec_split[2]) #Get just the species name
  cat("Calculating polygon areas for",spec_name,".\n")
  data <- read.csv(species) #Load the occurrence data back in
  pts <- terra::vect(data, geom = c("Longitude","Latitude"), crs = "epsg:7844") #Convert into a SpatVect object - spatial points
  #Also sets the Coordinate Reference System (EPSG:7844) - this is a suitable projection for Australian data
  
  #This approach creates a "buffer" around the points to estimate size
  #It is worth noting that this is /not/ the recommended strategy from the Threatened Species Scientific Committee
  buff <- as.polygons(terra::buffer(pts, width = buffer_width)) #Create a buffer of 2km (20000 m) around each point
  simp <- terra::aggregate(buff)
  area_buff <- round(terra::expanse(simp, unit = "km"), 2)
  
  #Minimum convex polygon approach
  #This is the standard way to calculate EOO for assessment
  mcp <- convHull(pts)
  #Clip to the coastline
  writeVector(mcp, filename = paste0(spec_name,"_minimum_convex_polygon.shp"), overwrite = TRUE) #Save the minimum convex polygon for future use
  area <- round(terra::expanse(mcp, unit = "km"), 2)
  cat("Extent of Occurrence = ",area,"km\U00B2.\n")
  
  #Area of Occupancy
  # Define a raster grid with 2x2 km cell size
  ext <- ext(simp)  # Get the extent of the merged buffer area
  raster_grid <- rast(ext, resolution = 0.021972184822683094,
                      crs = "EPSG:7844") #This is 2 x 2 km at the relevant latitude
  #Rasterize the points (cells will have a value of 1 where points exist, otherwise NA)
  point_raster <- rasterize(pts, raster_grid, fun = "count")
  
  #Determine the number of cells that contain at least one point
  occupied_cells <- sum(!is.na(values(point_raster)))
  
  #Calculate the Area of Occupancy (AOO) in square kilometers (2x2 km cells = 4 km² per cell)
  aoo_km2 <- occupied_cells * 4
  cat("Area of Occupancy = ",aoo_km2,"km\U00B2.\n")
  
  spec_df <- data.frame(spec_name, area, aoo_km2) #Combine area results together
  poly_df <- as.data.frame(rbind(poly_df, spec_df))  #Append the area results into the previous object
}

colnames(poly_df) <- c("Species","EOO","AOO") #Rename columns
write.csv(x = poly_df, file = "Marianthus_combined_area_estimates.csv", quote = F, row.names = F)

####Calculate the overlap of all pairs of polygons####
polygons <- list.files(path = "./", pattern = ".shp", full.names = F) #Get a list of the EOO polygons

#Initialize an empty SpatVector to store the combined polygons
combined_vector <- vect()

# Loop through the list of shapefiles and add each polygon shapefile
for (shp in polygons) {
  temp_vector <- vect(shp)  # Load the shapefile
  # Add the filename as an attribute/column in the SpatVector
  temp_vector$source_file <- basename(shp)  # Just the filename (without path)
  
  #Some of the polygons are empty due to too few occurrences - skip these ones as they break the code
  # heck if the geometry type is "polygons" before combining
  if (geomtype(temp_vector) == "polygons") {
    combined_vector <- rbind(combined_vector, temp_vector)
  } else {
   cat("Skipping", shp, "due to incompatible geometry type.\n")
  }
}

#Initialize a matrix to store the overlap areas
n <- length(combined_vector$source_file)  # Number of polygons
overlap_matrix <- matrix(0, n, n, dimnames = list(1:n, 1:n))
row.names(overlap_matrix) <- combined_vector$source_file #Set row names to the file names
colnames(overlap_matrix) <- combined_vector$source_file #Set column names to the file names

# Loop through each pair of polygons
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    # Calculate the intersection of the two polygons
    intersection <- terra::intersect(combined_vector[i, ], combined_vector[j, ])
    
    #If there is an intersection, calculate the area
    if (!is.empty(intersection)) {
      overlap_area <- expanse(intersection, unit = "km")  #Area in square kilometers
      overlap_matrix[i, j] <- overlap_area
      overlap_matrix[j, i] <- overlap_area  #Symmetric matrix
    }
  }
}

#View the overlap matrix (areas in square kilometers)
print(overlap_matrix)

#Convert to a paired list for simplicity
overlap_list <- melt(overlap_matrix)
#Simplify the file names to species
overlap_list$Var1 <- gsub(pattern = "_minimum_convex_polygon.shp", replacement = "", x = overlap_list$Var1)
overlap_list$Var2 <- gsub(pattern = "_minimum_convex_polygon.shp", replacement = "", x = overlap_list$Var2)
colnames(overlap_list) <- c("Species1", "Species2","OverlapArea_km2")

write.csv(x = overlap_list, file = "Marianthus_overlap_area_km2.csv", row.names = F, quote = F)

####Calculate percentage of range within overlap####
comb <- overlap_list %>%
  left_join(y = poly_df, by = join_by("Species1"=="Species")) %>%
  left_join(y = poly_df, by = join_by("Species2"=="Species")) %>%
  mutate(overlap.perc.s1 = round((OverlapArea_km2/EOO.x)*100, 2),
         overlap.perc.s2 = round((OverlapArea_km2/EOO.y)*100, 2))

colnames(comb) <- c("Species1","Species2","Overlap_Area_km2","Species1_EOO","Species1_AOO","Species2_EOO","Species2_AOO","Species1_Perc_EOO_Overlap","Species2_EOO_Range_Overlap")
write.csv(x = comb, file = "Marianthus_overlap_area_percentages.csv", row.names = FALSE, quote = F)
