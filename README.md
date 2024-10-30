# Marianthus hybridisation student project
Repository code for student project on hybridisation risks in Marianthus.

### Marianthus_ALA_analysis.R 
This R code downloads publicly available occurrence data from the Atlas of Living Australia (https://www.ala.org.au/) for all described species of the flowering plant genus Marianthus, and conducts several spatial analyses. These include:
1) Removing spatial outliers, which are >3 standard deviations from the mean latitude and longitude (to correct anomalies).
2) For each species, estimating the Extent of Occurrence (EOO) using a minimum convex polygon, following IUCN guidelines (https://www.ala.org.au/spatial-portal-help/aoo/)
3) For each species, estimating the Area of Occupancy (AOO) using a 2km x 2km raster, following IUCN guidelines (https://www.ala.org.au/spatial-portal-help/aoo/)
4) For each pair of species, calculating the absolute (km2) and relative (% of range) overlap in distribution based on the minimum convex polygons used to calculate EOO.

These calculations are output as multiple CSV files containing the relevant data.

