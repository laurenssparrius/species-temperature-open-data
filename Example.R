##########################################################################################
# STI calculation 
# v 1.0 Gregory van der Top - original version
# v 1.1 Laurens Sparrius - translated, csv support, clean up
# v 1.2 Dion van der Hak - added percentile calculations, improved performance, improved readability
# v 1.3 Dion van der Hak - rewrite: greatly improved performance by vectorization
#
# required files:
# 1. GBIF csv: e.g. gbif_plants.CSV (fields: species, decimallongitude, decimallatitude)
# 2. 50 km UTM grid: Grid_LAEA5210_50K_polygons.shp 
#    source: http://www.eea.europa.eu/data-and-maps/data/eea-reference-grids
# 3. BioClim dataset: C:\\avgtemppergridcelleurope
##########################################################################################

rm(list = ls()); gc()

#load functions
source("STI_plant_1_50km_per_250km.R")

occurrenceData <- read.csv("Xanthoria (test input).csv")
coordinateSystem <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

library(sp)
library(rgdal)
#read shapefile 50x50km grid Europe
#this shapefile has a Mercator projection as CRS
raster <- readOGR(dsn="Grid_ETRS89_LAEA5210_50KEEA15975I",layer="Grid_LAEA5210_50K_polygons")
meanTemperature <- readRDS("avgtemppergridcelleurope")

#execute the function
STISxx <- STISfunction(occurrenceData = occurrenceData, raster = raster, meanTemperature = meanTemperature, coordinateSystem = coordinateSystem)

head(STISxx)
write.csv(STISxx, file = "Xanthoria STI (test output).csv")

###############################