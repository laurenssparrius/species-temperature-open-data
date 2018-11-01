
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

STISfunction <- function(occurrenceData, raster, meanTemperature, coordinateSystem, nLoop = 100) {
  
  #stop function if wrong arguments are supplied
  stopifnot(is.numeric(occurrenceData$decimallatitude))
  stopifnot(is.numeric(occurrenceData$decimallongitude))
  stopifnot(is.factor(occurrenceData$species))
  stopifnot(is.factor(meanTemperature$CellCode))
  stopifnot(is.numeric(meanTemperature$EofOrigin))
  stopifnot(is.numeric(meanTemperature$NofOrigin))
  stopifnot(is.numeric(meanTemperature[,4]))
  stopifnot(is.character(coordinateSystem))
  stopifnot(is.numeric(nLoop))
  
  #load packages
  require(sp)
  require(data.table)
  require(rgdal)
  
  occurrenceData = data.table(occurrenceData)
  
  #register coordinate fields and coordinate reference system
  coordinates(occurrenceData) <- ~decimallongitude+decimallatitude
  proj4string(occurrenceData) <- CRS(coordinateSystem)
  
  #transform occurrence point data to Mercator projection
  occurrences <- spTransform(occurrenceData, proj4string(raster))
  occurrencesCoords <- occurrences@data
  occurrencesAttrib <- occurrences@coords
  occurrences <- data.table(occurrencesCoords, occurrencesAttrib)
  
  #round coordinates to match the 50 km grid
  occurrences$xhokcor <- mfloor((occurrences$decimallongitude / 10000), 5)
  occurrences$yhokcor <- mfloor((occurrences$decimallatitude / 10000), 5)
  
  #associate grid cells with cell codes of the ETRS grid map
  occurrences$Cellcode <- paste0("50kmE", occurrences$xhokcor, "N", occurrences$yhokcor)
  
  #add generalized 250 km grid cells on top
  occurrences$X250kmhok <- mfloor(occurrences$xhokcor, 25)
  occurrences$y250kmhok <- mfloor(occurrences$yhokcor, 25)
  occurrences$Cell250 <- paste0("250kmE", occurrences$X250kmhok, "N", occurrences$y250kmhok)
  
  #remove unused data
  occurrences <- occurrences[, .(species, Cellcode, Cell250)]
  meanTemperature <- data.table(meanTemperature[meanTemperature$CellCode %in% occurrences$Cellcode,])
  
  #50 - 250 cell associations
  cellcodes <- occurrences[!duplicated(occurrences$Cellcode)]
  cellcodes <- cellcodes[, .(Cellcode, Cell250)]
  #temperature
  cellcodes <- merge(cellcodes, meanTemperature[, c(1,4)], by.x = "Cellcode", by.y = "CellCode")
  cellcodes$gemtempeuropa.1 <- cellcodes$gemtempeuropa.1 / 10
  cellcodes$Cellcode <- as.character(cellcodes$Cellcode)
  colnames(cellcodes) <- c("Cellcode", "Cell250", "temperature")
  
  arrBootstrapSTI <- NULL
  
  for(i in 1:nLoop) {
    #random select 1 50km cell within each 250km cell
    randomcells <- tapply(cellcodes$Cellcode, cellcodes$Cell250, sample, size = 1)
    
    ##############################
    #Calculate STI for all species
    ##############################
    
    #retrieve occurrences for each species
    plant250 <- subset(occurrences, Cellcode %in% randomcells)
    #add temperature per occurrence
    STIdata <- merge(plant250, cellcodes[, c(1,3)], by = "Cellcode")
    
    #Calculate STI from array
    STI <- tapply(STIdata$temperature, STIdata$species, FUN = mean, na.rm = T)
    STI <- data.frame(STI)
    STI <- data.table(species = row.names(STI), STI = STI$STI)
    
    stdv <- tapply(STIdata$temperature, STIdata$species, FUN = sd, na.rm = T)
    stdv <- data.frame(stdv)
    stdv <- data.table(species = row.names(stdv), stdv = stdv$stdv)
    STIdatawithspecies <- merge(STI, stdv, by = "species")
    
    perc5 <- tapply(STIdata$temperature, STIdata$species, FUN = quantile, probs = 0.05, na.rm = T)
    perc5 <- data.frame(perc5)
    perc5 <- data.table(species = row.names(perc5), perc5 = perc5$perc5)
    STIdatawithspecies <- merge(STIdatawithspecies, perc5, by = "species")
    
    perc25 <- tapply(STIdata$temperature, STIdata$species, FUN = quantile, probs = 0.25, na.rm = T)
    perc25 <- data.frame(perc25)
    perc25 <- data.table(species = row.names(perc25), perc25 = perc25$perc25)
    STIdatawithspecies <- merge(STIdatawithspecies, perc25, by = "species")
    
    perc75 <- tapply(STIdata$temperature, STIdata$species, FUN = quantile, probs = 0.75, na.rm = T)
    perc75 <- data.frame(perc75)
    perc75 <- data.table(species = row.names(perc75), perc75 = perc75$perc75)
    STIdatawithspecies <- merge(STIdatawithspecies, perc75, by = "species")
    
    perc95 <- tapply(STIdata$temperature, STIdata$species, FUN = quantile, probs = 0.95, na.rm = T)
    perc95 <- data.frame(perc95)
    perc95 <- data.table(species = row.names(perc95), perc95 = perc95$perc95)
    STIdatawithspecies <- merge(STIdatawithspecies, perc95, by = "species")
    
    arrBootstrapSTI <- rbind(arrBootstrapSTI, STIdatawithspecies)
  }
  
  #get the averages over the samples
  arrFinalSTI <- tapply(arrBootstrapSTI$STI, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalSTI <- data.frame(arrFinalSTI)
  arrFinalSTI <- data.table(species = row.names(arrFinalSTI), STI = arrFinalSTI$arrFinalSTI)
  
  arrFinalStDev <- tapply(arrBootstrapSTI$stdv, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalStDev <- data.frame(arrFinalStDev)
  arrFinalStDev <- data.table(species = row.names(arrFinalStDev), SD = arrFinalStDev$arrFinalStDev)
  STISxx <- merge(arrFinalSTI, arrFinalStDev, by = "species")
  
  arrFinalPerc5 <- tapply(arrBootstrapSTI$perc5, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalPerc5 <- data.frame(arrFinalPerc5)
  arrFinalPerc5 <- data.table(species = row.names(arrFinalPerc5), perc5 = arrFinalPerc5$arrFinalPerc5)
  STISxx <- merge(STISxx, arrFinalPerc5, by = "species")
  
  arrFinalPerc25 <- tapply(arrBootstrapSTI$perc25, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalPerc25 <- data.frame(arrFinalPerc25)
  arrFinalPerc25 <- data.table(species = row.names(arrFinalPerc25), perc25 = arrFinalPerc25$arrFinalPerc25)
  STISxx <- merge(STISxx, arrFinalPerc25, by = "species")
  
  arrFinalPerc75 <- tapply(arrBootstrapSTI$perc75, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalPerc75 <- data.frame(arrFinalPerc75)
  arrFinalPerc75 <- data.table(species = row.names(arrFinalPerc75), perc75 = arrFinalPerc75$arrFinalPerc75)
  STISxx <- merge(STISxx, arrFinalPerc75, by = "species")
  
  arrFinalPerc95 <- tapply(arrBootstrapSTI$perc95, arrBootstrapSTI$species, FUN = mean, na.rm = T)
  arrFinalPerc95 <- data.frame(arrFinalPerc95)
  arrFinalPerc95 <- data.table(species = row.names(arrFinalPerc95), perc95 = arrFinalPerc95$arrFinalPerc95)
  STISxx <- merge(STISxx, arrFinalPerc95, by = "species")
  
  #compile dataset for export
  colnames(STISxx)[4:7] <- c("5th_percentile", "25th_percentile", "75th_percentile", "95th_percentile")
  
  return(STISxx)
}

#function to round to the nearest base value
mfloor <- function(x, base) {
  base * floor(x / base)
}
