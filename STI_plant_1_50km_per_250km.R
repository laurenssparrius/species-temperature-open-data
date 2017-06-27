##########################################################################################
# STI calculation 
# v 1.0 Gregory van der Top - original version
# v 1.1 Laurens Sparrius - translated, csv support, clean up
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

arrSti<-NULL

occurrenceData<-read.csv("C:\\Users\\Laurens\\Desktop\\plantenvoorpub\\gbif_data.csv")  
data.table(occurrenceData)

#register coordinate fields and coordinate reference system
coordinates(occurrenceData)<-~decimallongitude+decimallatitude
proj4string(occurrenceData)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#read shapefile 50x50km grid Europe
#this shapefile has a Mercator projection as CRS
rasterEurope<-readOGR(dsn="C:/Users/Laurens/Desktop/plantenvoorpub/Grid_ETRS89_LAEA5210_50KEEA15975I",layer="Grid_LAEA5210_50K_polygons")

#transform occurrence point data to Mercator projection
occurrences<-spTransform(occurrenceData,proj4string(rasterEurope))
occurrencesCoords<-occurrences@data
occurrencesAttrib<-occurrences@coords
occurrences<-data.table(occurrencesCoords,occurrencesAttrib)

#function to round to the nearest base value
mfloor<-function(x,base){
  base*floor(x/base)
}

#round coordinates to match the 50 km grid
occurrences$xhokcor<-mfloor((occurrences$decimallongitude/10000),5)
occurrences$yhokcor<-mfloor((occurrences$decimallatitude/10000),5)

#associate grid cells with cell codes of the ETRS grid map
occurrences$"50kmE"<-"50kmE"
occurrences$N<-"N"
occurrences$Cellcode<-paste(occurrences$`50kmE`,occurrences$xhokcor,occurrences$N,occurrences$yhokcor, sep = "")
occurrences$CellCode<-factor(occurrences$Cellcode)


#add generalized 250 km grid cells on top
occurrences$X250kmhok<-mfloor(occurrences$xhokcor, 25)
occurrences$y250kmhok<-mfloor(occurrences$yhokcor, 25)
occurrences$"250kmE"<-"250kmE"
occurrences$Cell250<-paste(occurrences$`250kmE`,occurrences$X250kmhok,occurrences$N,occurrences$y250kmhok, sep = "")

  #open climate data
  meanTemperatureEurope<-readRDS("C:\\Users\\Laurens\\Desktop\\plantenvoorpub\\avgtemppergridcelleurope")

STISfunction<-function(x){

  arrStiStd<-NULL

  #random select 1 50km cell within each 250km cell
  for (i in unique(x$Cell250)){
    gridcell250<-x[x$Cell250==i,]
    randomgridcell<-sample(unique(gridcell250$CellCode), 1)
    plant250<-subset(x, CellCode %in% randomgridcell)
    arrStiStd<-rbind(arrStiStd, plant250)
   }
  
  plant250<-arrStiStd
  
  ############################
  #Calculate STI for a species
  ############################

  arrSti<-NULL
  arrStiStd<-NULL

  for (i in unique(plant250$species)){

    #retrieve occurrences for each species
    sampleOccurrences<-plant250[species==i]

    STIdata<-NULL

    for (i in unique(sampleOccurrences$CellCode)){
      tempx<-meanTemperatureEurope$gemtempeuropa.1[meanTemperatureEurope$CellCode==i]
      STIdata<-rbind(STIdata,tempx)
      STIdata2<-STIdata/10
    }
    
    #Calculate STI from array
    STI<-mean(STIdata2, na.rm=T)
    stdv<-sd(STIdata2, na.rm=T)
    arrSti<-rbind(arrSti, STI)
    arrStiStd<-rbind(arrStiStd, stdv)
    
  }
  
  colnames(arrSti)[1]<-"STIspecies"
  colnames(arrStiStd)[1]<-"STDEVspecies"
  STIdata<-data.table(arrSti)
  STI.STDEV<-data.table(arrStiStd)

  speciesname<-unique(plant250$species)
  STIdatawithspecies<-data.table(STIdata, speciesname, STI.STDEV)
  
}

arrBootstrapSti<-NULL

for(i in 1:100){
  STIS<-STISfunction(x=occurrences)
  arrBootstrapSti<-rbind(arrBootstrapSti, STIS)
}

arrFinalSTI<-NULL
arrFinalStDev<-NULL

for(i in unique(arrBootstrapSti$speciesname)){
  uniqueSpecies<-arrBootstrapSti[arrBootstrapSti$speciesname==i,]
  speciesname=i
  STImean<-mean(uniqueSpecies$STIspecies, na.rm=T)
  STIsddev<-mean(uniqueSpecies$STDEVspecies,na.rm=T)
  arrFinalSTI<-rbind(arrFinalSTI, STImean)
  arrFinalStDev<-rbind(meanstdev,STIsddev)
}

#compile dataset for export
STISmean<-data.table(arrFinalSTI)
STISstdev<-data.table(arrFinalStDev)
colnames(STISmean)[1]<-"STIspecies"
colnames(STISstdev)[1]<-"STIstdev"
STISxx<-data.table(STISmean,STISstdev)
STISxx$species<-unique(arrBootstrapSti$speciesname)
colnames(STISxx)[1]<-"STI"
colnames(STISxx)[2]<-"SD"

write.csv(STISxx, file = "C:\\Users\\Laurens\\Desktop\\plantenvoorpub\\SpeciesSTI.csv")
##############################################
