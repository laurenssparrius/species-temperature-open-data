##########################################################################################
# STI calculation 
# v 1.0 Gregory van der Top - original version
# v 1.1 Laurens Sparrius - translated, csv support, clean up
# v 1.2 Dion van der Hak - added percentile calculations, improved performance, improved readability, added multithreading support
#
# required files:
# 1. GBIF csv: e.g. gbif_plants.CSV (fields: species, decimallongitude, decimallatitude)
# 2. 50 km UTM grid: Grid_LAEA5210_50K_polygons.shp 
#    source: http://www.eea.europa.eu/data-and-maps/data/eea-reference-grids
# 3. BioClim dataset: C:\\avgtemppergridcelleurope
##########################################################################################

rm(list=ls())
library(sp)
library(data.table)
library(rgdal)

arrSti <- NULL

occurrenceData <- read.csv("\\species-temperature-open-data\\Xanthoria (test input).csv")  
data.table(occurrenceData)

#register coordinate fields and coordinate reference system
coordinates(occurrenceData) <- ~decimallongitude+decimallatitude
proj4string(occurrenceData) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#read shapefile 50x50km grid Europe
#this shapefile has a Mercator projection as CRS
rasterEurope <- readOGR(dsn="\\species-temperature-open-data\\Grid_ETRS89_LAEA5210_50KEEA15975I",layer="Grid_LAEA5210_50K_polygons")

#transform occurrence point data to Mercator projection
occurrences <- spTransform(occurrenceData,proj4string(rasterEurope))
occurrencesCoords <- occurrences@data
occurrencesAttrib <- occurrences@coords
occurrences <- data.table(occurrencesCoords,occurrencesAttrib)

#function to round to the nearest base value
mfloor <- function(x, base) {
  base * floor(x / base)
}

#round coordinates to match the 50 km grid
occurrences$xhokcor <- mfloor((occurrences$decimallongitude / 10000), 5)
occurrences$yhokcor <- mfloor((occurrences$decimallatitude / 10000), 5)

#associate grid cells with cell codes of the ETRS grid map
occurrences$"50kmE" <- "50kmE"
occurrences$N <- "N"
occurrences$Cellcode <- paste0(occurrences$`50kmE`, occurrences$xhokcor, occurrences$N, occurrences$yhokcor)

#add generalized 250 km grid cells on top
occurrences$X250kmhok <- mfloor(occurrences$xhokcor, 25)
occurrences$y250kmhok <- mfloor(occurrences$yhokcor, 25)
occurrences$"250kmE" <- "250kmE"
occurrences$Cell250 <- paste0(occurrences$`250kmE`, occurrences$X250kmhok, occurrences$N, occurrences$y250kmhok)

#open climate data
meanTemperatureEurope <- readRDS("\\species-temperature-open-data\\avgtemppergridcelleurope")
meanTemperatureEurope$CellCode <- as.character(meanTemperatureEurope$CellCode)
meanTemperatureEurope <- meanTemperatureEurope[meanTemperatureEurope$CellCode %in% occurrences$Cellcode,] #remove unused temperature data to speed up calculations

STISfunction <- function(x) {
  
  x$species <- as.character(x$species)

  arrStiStd <- NULL

  #random select 1 50km cell within each 250km cell
  for(i in unique(x$Cell250)) {
    gridcell250 <- x[x$Cell250 == i,]
    randomgridcell <- sample(unique(gridcell250$Cellcode), 1)
    plant250 <- subset(x, Cellcode %in% randomgridcell)
    arrStiStd <- rbind(arrStiStd, plant250)
   }
  
  plant250 <- arrStiStd
  
  ############################
  #Calculate STI for a species
  ############################

  arrSti <- NULL
  arrStiStd <- NULL
  arrStiPerc5 <- NULL
  arrStiPerc25 <- NULL
  arrStiPerc75 <- NULL
  arrStiPerc95 <- NULL

  for(i in unique(plant250$species)) {

    #retrieve occurrences for each species
    sampleOccurrences <- plant250[species == i]

    STIdata <- NULL

    for(j in unique(sampleOccurrences$Cellcode)) {
      tempx <- meanTemperatureEurope$gemtempeuropa.1[meanTemperatureEurope$CellCode == j]
      STIdata <- rbind(STIdata, tempx)
      STIdata2 <- STIdata / 10
    }
    
    #Calculate STI from array
    STI <- mean(STIdata2, na.rm = T)
    stdv <- sd(STIdata2, na.rm = T)
    perc5 <- quantile(STIdata2, 0.05, na.rm = T)
    perc25 <- quantile(STIdata2, 0.25, na.rm = T)
    perc75 <- quantile(STIdata2, 0.75, na.rm = T)
    perc95 <- quantile(STIdata2, 0.95, na.rm = T)
    
    arrSti <- rbind(arrSti, STI)
    arrStiStd <- rbind(arrStiStd, stdv)
    arrStiPerc5 <- rbind(arrStiPerc5, perc5)
    arrStiPerc25 <- rbind(arrStiPerc25, perc25)
    arrStiPerc75 <- rbind(arrStiPerc75, perc75)
    arrStiPerc95 <- rbind(arrStiPerc95, perc95)
  }
  
  colnames(arrSti)[1] <- "STIspecies"
  colnames(arrStiStd)[1] <- "STDEVspecies"
  colnames(arrStiPerc5)[1] <- "perc5species"
  colnames(arrStiPerc25)[1] <- "perc25species"
  colnames(arrStiPerc75)[1] <- "perc75species"
  colnames(arrStiPerc95)[1] <- "perc95species"
  STIdata <- data.table(arrSti)
  STI.STDEV <- data.table(arrStiStd)
  STI.perc5 <- data.table(arrStiPerc5)
  STI.perc25 <- data.table(arrStiPerc25)
  STI.perc75 <- data.table(arrStiPerc75)
  STI.perc95 <- data.table(arrStiPerc95)

  speciesname <- unique(plant250$species)
  STIdatawithspecies <- data.table(STIdata, speciesname, STI.STDEV, STI.perc5, STI.perc25, STI.perc75, STI.perc95)
  
  return(STIdatawithspecies)
}

#multithreaded stuff
library(foreach)
library(doParallel)
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

arrBootstrapSti <- NULL
nLoop <- 100
time1 <- as.numeric(Sys.time())

foreach(i = 1:nLoop, .combine = rbind, .packages = 'data.table') %dopar% {
  STIS <- STISfunction(x = occurrences)
  STIS
}

time2 <- as.numeric(Sys.time())
print(paste0("Finished in ", time2 - time1, " seconds!"))

stopCluster(cl)

arrBootstrapSti <- STIS

arrFinalSTI <- NULL
arrFinalStDev <- NULL
arrFinalPerc5 <- NULL
arrFinalPerc25 <- NULL
arrFinalPerc75 <- NULL
arrFinalPerc95 <- NULL

for(i in unique(arrBootstrapSti$speciesname)) {
  uniqueSpecies <- arrBootstrapSti[arrBootstrapSti$speciesname == i,]
  speciesname <- i
  STImean <- mean(uniqueSpecies$STIspecies, na.rm = T)
  STIsddev <- mean(uniqueSpecies$STDEVspecies, na.rm = T)
  STIperc5 <- mean(uniqueSpecies$perc5species, na.rm = T)
  STIperc25 <- mean(uniqueSpecies$perc25species, na.rm = T)
  STIperc75 <- mean(uniqueSpecies$perc75species, na.rm = T)
  STIperc95 <- mean(uniqueSpecies$perc95species, na.rm = T)
  arrFinalSTI <- rbind(arrFinalSTI, STImean)
  arrFinalStDev <- rbind(arrFinalStDev, STIsddev)
  arrFinalPerc5 <- rbind(arrFinalPerc5, STIperc5)
  arrFinalPerc25 <- rbind(arrFinalPerc25, STIperc25)
  arrFinalPerc75 <- rbind(arrFinalPerc75, STIperc75)
  arrFinalPerc95 <- rbind(arrFinalPerc95, STIperc95)
}

#compile dataset for export
STISmean <- data.table(arrFinalSTI)
STISstdev <- data.table(arrFinalStDev)
STISperc5 <- data.table(arrFinalPerc5)
STISperc25 <- data.table(arrFinalPerc25)
STISperc75 <- data.table(arrFinalPerc75)
STISperc95 <- data.table(arrFinalPerc95)
colnames(STISmean)[1] <- "STIspecies"
colnames(STISstdev)[1] <- "STIstdev"
colnames(STISperc5)[1] <- "5th_percentile"
colnames(STISperc25)[1] <- "25th_percentile"
colnames(STISperc75)[1] <- "75th_percentile"
colnames(STISperc95)[1] <- "95th_percentile"
STISxx <- data.table(STISmean, STISstdev, STISperc5, STISperc25, STISperc75, STISperc95)
STISxx$species <- unique(arrBootstrapSti$speciesname)
colnames(STISxx)[1] <- "STI"
colnames(STISxx)[2] <- "SD"

write.csv(STISxx, file = "\\species-temperature-open-data\\SpeciesSTI.csv")
##############################################
