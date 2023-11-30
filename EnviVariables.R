
library(terra)
library(sf)
library(lubridate)
library(sp)
library(magrittr)
library(ggplot2)
library(dplyr)  
library(reshape2)
## -------load shapefiles-------------


# tracks for counting seabirds between Punta Arenas and King George Island
#30 to 40 min transects registering all birds within 300 to 500m from ship 
#1h minimal between transects 
#birds following ship were counted once when entering the transect


# during transects a hand held GPS was set to auto-tracking, marking one point each 5 seconds. 

t20<-st_as_sf(st_read("C:/OnboardSeabirdCounts/Track_ECA57.shp"))
t21<-st_as_sf(st_read("C:/OnboardSeabirdCounts/Track_ECA58.shp"))
t22<-st_as_sf(st_read("C:/OnboardSeabirdCounts/Track_ECA59b.shp"))

head(t22)

### geographical coordinates for each geographical fix


coordf20 <- t20 %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

coordf21<- t21 %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

coordf22 <- t22 %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])


# time stamps 

coordf20$timestamp<-as.POSIXct(strptime(coordf20$DateTimeS, format="%Y-%m-%d T %H:%M:%S", tz="GMT"))
coordf21$timestamp<-as.POSIXct(strptime(coordf21$DateTimeS, format="%Y-%m-%d T %H:%M:%S", tz="GMT"))
coordf22$timestamp<-as.POSIXct(strptime(coordf22$BeginTime, format="%Y-%m-%d %H:%M:%S", tz="GMT")) # this year I used a Garmin watch, so the configuration was different

summary(coordf20$timestamp)

# standardize data sets
head(coordf20)

cdf20<-data.frame(coordf20[12],coordf20[14:16])
cdf21<-data.frame(coordf21[12],coordf21[14:16])

head(coordf22)
cdf22<-data.frame(coordf22[15],coordf22[17:19])


transects<-rbind(cdf20,cdf21,cdf22)

transects$Year<-year(transects$timestamp)
transects$Month<-month(transects$timestamp)
transects$Day<-day(transects$timestamp)
transects$Hour<-hour(transects$timestamp)

head(transects)

trs<-transects

head(trs)

summary(as.factor(trs$Year))


trsm<-plyr::ddply(trs, c("Transect","Year","Month","Day","Hour"), summarise,
                  lon.min=min(lon),lon.max=max(lon),
                  lat.min=min(lat),max.lat=max(lat),
                  lon=mean(lon),lat=mean(lat))

trsm$timestamp<-as.POSIXct(strptime(paste(paste(trsm$Year,trsm$Month,trsm$Day,sep="-"),
                                          trsm$Hour), format="%Y-%m-%d %H", tz="GMT"))


### download environmental data from Copernicus Marine Services

#https://data.marine.copernicus.eu/viewer



library(matrixStats)

library(reshape2)
library(ggplot2)
library(lubridate)

library(raster)
library(terra)
library(sf)
library(dplyr)
library(ncdf4)

library(raster)
library(ncdf4)

library(raster)
library(ncdf4)


### extract variables into transects based on the timestamp

###------------Wind variables-----------------

# product ID WIND_GLO_PHY_CLIMATE_L4_MY_012_003
#https://doi.org/10.48670/moi-00181


#### ------eastward wind (per month) ---------------
ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_obs-wind_glo_phy_my_l4_P1M_1701281468892"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "eastward_wind"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("Ewind_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "hours")) + 1)/12
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_Ewind <- do.call(rbind, extracted_values)

head(result_Ewind)
result_Ewind<-data.frame(result_Ewind)

trsm$Ewind<-result_Ewind$value


ggplot(trsm,aes(lat,Ewind,colour=as.factor(Year)))+geom_point()


#### ------northward wind (per month) ---------------
ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_obs-wind_glo_phy_my_l4_P1M_1701281468892"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "northward_wind"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]

    # Save raster to file
  writeRaster(tmp_raster, filename = paste0("Nwind_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)/12
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_Nwind <- do.call(rbind, extracted_values)

head(result_Nwind)
result_Nwind<-data.frame(result_Nwind)

trsm$Nwind<-result_Nwind$value

trsm$windspeed<-sqrt((trsm$Ewind*trsm$Ewind)+(trsm$Nwind*trsm$Nwind))

ggplot(trsm,aes(lat,windspeed,colour=as.factor(Year)))+geom_point()

####------------CHL per day-------------

#Product ID: GLOBAL_ANALYSISFORECAST_BGC_001_028
#https://doi.org/10.48670/moi-00015

ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "global-analysis-forecast-bio-001-028-daily_1701280896260"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "chl"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("chl_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}


# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_chl <- do.call(rbind, extracted_values)

head(result_chl)
result_chl<-data.frame(result_chl)

trsm$chl<-result_chl$value


ggplot(trsm,aes(lat,chl,colour=as.factor(Year)))+geom_point()


###------- SST -------------also per day

#Product ID : GLOBAL_ANALYSISFORECAST_PHY_001_024
#https://doi.org/10.48670/moi-00016

ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m_1701280484714"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "thetao"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("sst_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_sst <- do.call(rbind, extracted_values)

head(result_sst)
result_sst<-data.frame(result_sst)

trsm$sst<-result_sst$value


ggplot(trsm,aes(lat,sst,colour=as.factor(Year)))+geom_point()





### -----------ocean mixed layer thickness defined by sigma------------
## product ID GLOBAL_ANALYSISFORECAST_PHY_001_024
#https://doi.org/10.48670/moi-00016

ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_mod_glo_phy_anfc_0.083deg_P1D-m_1701280574463"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "mlotst"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("mld_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_mld <- do.call(rbind, extracted_values)

head(result_mld)
result_mld<-data.frame(result_mld)

trsm$mld<-result_mld$value


ggplot(trsm,aes(lat,mld,colour=as.factor(Year)))+geom_point()


###---------- sea water velocity eastward --------------------

# product ID: GLOBAL_ANALYSISFORECAST_PHY_001_024

#https://doi.org/10.48670/moi-00016


ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m_1701280539602"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "uo"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("uo_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_uo<- do.call(rbind, extracted_values)

head(result_uo)
result_uo<-data.frame(result_uo)

trsm$uo<-result_uo$value


ggplot(trsm,aes(lat,uo,colour=as.factor(Year)))+geom_point()




###---------- sea water velocity northward --------------------

# product ID: GLOBAL_ANALYSISFORECAST_PHY_001_024

#https://doi.org/10.48670/moi-00016


ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m_1701280539602"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "vo"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("vo_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_vo<- do.call(rbind, extracted_values)

head(result_vo)
result_vo<-data.frame(result_vo)

trsm$vo<-result_vo$value


ggplot(trsm,aes(lat,vo,colour=as.factor(Year)))+geom_point()


trsm$swv<-sqrt((trsm$uo*trsm$uo)+(trsm$vo*trsm$vo))

write.csv(trsm,"transects_envi_data.csv")

ggplot(trsm,aes(lat,swv,colour=as.factor(Year)))+geom_point()
ggplot(trsm,aes(sst,swv,colour=as.factor(Year)))+geom_point()


###---------- Mass content of zooplankton expressed as carbon in sea water (g/m2)--------------------

# product ID: GLOBAL_MULTIYEAR_BGC_001_033

#https://doi.org/10.48670/moi-00020


ncpath <- "C:/OnboardSeabirdCounts/EnviData/"
ncname <- "cmems_mod_glo_bgc_my_0.083deg-lmtl_PT1D-i_1701281199592"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "zooc"
tmp_raster <- brick(ncfname, varname=dname)
tmp_brick <- brick(ncfname, varname=dname)
plot(tmp_raster)


# Loop through each layer (day) in the brick
for (i in 1:nlayers(tmp_brick)) {
  # Extract the raster for the current day
  tmp_raster <- tmp_brick[[i]]
  
  
  # Save raster to file
  writeRaster(tmp_raster, filename = paste0("zooc_output_day_", i, ".tif"), format = "GTiff", overwrite = TRUE)
}



# Create an empty list to store the extracted values
extracted_values <- list()
tmp_brick
# Loop through each row in "trs"
for (i in 1:nrow(trsm)) {
  # Extract the timestamp, lat, and lon for the current row
  current_timestamp <- trsm$timestamp[i]
  current_lat <- trsm$lat[i]
  current_lon <- trsm$lon[i]
  
  # Find the corresponding raster layer in the brick based on the timestamp
  current_day <- as.integer(as.numeric(difftime(current_timestamp, min(trsm$timestamp), units = "days")) + 1)
  current_raster <- tmp_brick[[3]]
  
  # Create a SpatialPoints object for extraction
  points <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2), proj4string = CRS(proj4string(current_raster)))
  
  # Extract the value at the specified coordinates
  extracted_value <- extract(current_raster, points)
  
  # Store the extracted value in the list
  extracted_values[[i]] <- c(timestamp = current_timestamp, value = extracted_value)
}

# Convert the list to a data frame
result_zooc<- do.call(rbind, extracted_values)

head(result_zooc)
result_zooc<-data.frame(result_zooc)

trsm$zooc<-result_zooc$value


ggplot(trsm,aes(lat,zooc,colour=as.factor(Year)))+geom_point()


write.csv(trsm,"transects_envi_data.csv")

ggplot(trsm,aes(lat,swv,colour=as.factor(Year)))+geom_point()
ggplot(trsm,aes(sst,swv,colour=as.factor(Year)))+geom_point()
