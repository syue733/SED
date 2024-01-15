install.packages("chron")
install.packages("ncdf4")
library(ncdf4) 
library(sf)
library(tmap)
library(chron)

###predicted
solar <- nc_open("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/era6.nc")
solar
lon <- ncvar_get(solar, "longitude")
lat <- ncvar_get(solar, "latitude")
time <- ncvar_get(solar, "time")
time

dim(time)
tunits <- ncatt_get(solar,"time","units")


tustr <- strsplit(tunits$value, " ") 
tdstr <- strsplit(unlist(tustr)[3], "-") 
tyear <- as.integer(unlist(tdstr)[1]) 
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])

chron(time/24, origin=c(tmonth, tday, tyear) )
ssrd_array <- ncvar_get(solar,"ssrd")
ssrd_slice <- ssrd_array[,,2]

lonlat <- as.matrix( (expand.grid(lon, lat)))
dim(lonlat)
ssrd_vec <- as.vector( ssrd_slice) 
length(ssrd_vec)
ssrd_df <- data.frame( cbind( lonlat,ssrd_vec  ))
colnames(ssrd_df) <- c("lon", "lat", "ssrd")
ssrd_df_value <- na.omit (ssrd_df)


ssrd_sf<- st_as_sf( ssrd_df_value, coords = c(  "lon", "lat")  )
st_crs(ssrd_sf) <- 4326 
ssrd_sf <- st_transform(ssrd_sf, 4326 )

ssrd_sf = st_transform(ssrd_sf, 4326)
indonesia <- st_read("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/idn_admbnda_adm0_bps_20200401.shp")
indonesia = st_transform(indonesia, st_crs(ssrd_sf))

coor = as.data.frame(st_coordinates(ssrd_sf))
View(ssrd_sf)

ssrd_sf$x = coor$X
ssrd_sf$y = coor$Y
ssrd_nogeom = st_drop_geometry(ssrd_sf) #get rid of geometry but keep all other attributes
ssrd_nogeom=na.omit(ssrd_nogeom)
install.packages("gstat")

library(gstat)
gs <- gstat(formula=ssrd~1, locations=~x+y, data=ssrd_nogeom, nmax=Inf, set=list(idp=5)) #data should be in data frame format
gs
st_bbox(indonesia)
library(raster)
 
raster_template <- raster(extent(95.01079, 141.01940, -11.00762, 6.07693),
                          resolution = 0.05,
                          crs = st_crs(indonesia)$wkt)
raster_template
idw <- interpolate(raster_template, gs,z = ssrd_nogeom$ssrd, debug.level=0)
names(idw)                   
plot(idw$x)
idw_mask <- mask(idw, indonesia)
plot(idw_mask$x)
names(idw_mask) = c( "predicted_ssrd" )

tmap_mode("view")

tm_shape(idw_mask$predicted_ssrd) + 
  tm_raster(col="predicted_ssrd", style = "quantile", n = 10, palette= "Reds", legend.show = TRUE)+ tm_scale_bar(position = c("right", "bottom"))

###RMSE
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

null <- RMSE(mean(ssrd_sf$ssrd), ssrd_sf$wind)
null 

n_idp = 20 
n_fold =10

rmse <- rep(NA, n_fold) 
set.seed(7713)
kf <- sample(1:n_fold, nrow(ssrd_nogeom), replace=TRUE)
va = data.frame( c(1:n_idp), NA)
colnames(va) =c("idp","rmse") 

for (j in 1:n_idp) 
{
  for (i in 1:n_fold) {
    test <- ssrd_nogeom[kf == 1, ]
    train <- ssrd_nogeom[kf != 1, ]
    gs <- gstat(formula=ssrd~1, locations=~x+y, data=train, nmax=Inf, set=list(idp=j))
    pre = predict(gs, test, debug.level=0 )
    rmse[i] <- RMSE(test$ssrd, pre$var1.pred)
  }
  va[j,2] = (mean(rmse) )
}

va[which(va$rmse==min(va)),]

library(ggplot2)
ggplot(va) +
  geom_point(aes(x = idp, y= rmse))+
  geom_hline(yintercept=min(va), linetype="dashed", color = "red")+
  theme_classic()

###kwh
radiation_to_power <- function(G, A=1, r=0.175, p=0.6, hours=1){
  kWh <- G * A * r * p * (hours/3600) / 1000
  return(kWh)
}
idw_mask_kwh=radiation_to_power (idw_mask$predicted_ssrd)
names(idw_mask_kwh)=c("ssrd_kwh")

tm_shape(idw_mask_kwh)+
  tm_raster(col="ssrd_kwh", style = "quantile", n=10, palette = "YlOrRd",legend.show = TRUE)+
  tm_scale_bar(position = c("right", "bottom"))

radiation_to_power <- function(G, A=1, r=0.175, p=0.6, hours=1){
  kWh <- G * A * r * p * (hours/3600) / 1000
  return(kWh)
}
idw_mask_kwh=radiation_to_power (idw_mask$predicted_ssrd)
names(idw_mask_kwh)=c("ssrd_kwh")

tm <-tm_shape(idw_mask_kwh)+
  tm_raster(col="ssrd_kwh", style = "quantile", n=10, palette = "YlOrRd",legend.show = FALSE)+
  tm_scale_bar(position = c("right", "bottom"))

tmap_save(tm, filename = "/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/map.tiff",
          width = 10, height =10, units = "in", dpi = 300)

###Elevation
library(terra)
elevation_raster <- rast("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/IDN_alt.vrt")

indonesia <- vect("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/idn_admbnda_adm0_bps_20200401.shp")
crs(r)
crs(indonesia)
indonesia <- project(indonesia, crs(r))

r_masked <- mask(r, indonesia)
plot(r_masked)

tm_shape(r_masked) +
  tm_raster(col = "IDN_alt", style = "quantile", n = 10, palette = "Blues", legend.show = TRUE) +
  tm_view(legend.position = c("right", "bottom"))     

tmap_save(tm, filename = "/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/map.tiff",
          width = 10, height =10, units = "in", dpi = 300)

###slope
library(terra)

r <- rast("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/IDN_alt.vrt")

indonesia <- vect("/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/idn_admbnda_adm0_bps_20200401.shp")

indonesia <- project(indonesia, crs(r))

r_masked <- mask(r, indonesia)

slope <- terrain(r_masked, "slope")
slope_legend_labels <- c("0-1", "1-2", "2-3", "3-5", "5-10", "10-50")

plot(r_masked)


tm_shape(r_masked) +
  tm_raster(col = "IDN_alt", style = "quantile", n = 10, palette = "Blues", legend.show = FALSE) +
  tm_shape(slope) +
  tm_raster(col = "slope", style = "cont", title = "Slope") +
  tm_view(legend.position = c("right", "bottom"))

tmap_save(tm, filename = "/Users/bb/Desktop/0093-spatial analysis of ED/Week_6_Practical/Week_6_data/map.tiff",
          width = 10, height =10, units = "in", dpi = 300)


#1 Net Present Value (NPV)======
rep(10,4) #create a vector by repeating 10 4 times
# output of the function above is: 10 10 10 10
seq( 1, 11, 2) #create a sequence of data start from 1 and end at 11. 2 is the increment of the sequence.
#outout will be: 1 3 5 7 9 11

calc_NPV <- function(annual_revenue=8090000000, i=0.05, lifetime_yrs=25, CAPEX, OPEX=0){
  costs_op <- rep(OPEX, lifetime_yrs) #operating cost
  revenue <- rep(annual_revenue, lifetime_yrs) 
  t <- seq(1, lifetime_yrs, 1) #output: 1, 2, 3, ...25
  
  NPV <- sum( (revenue - costs_op)/(1 + i)**t ) - CAPEX
  return(round(NPV, 0))
}
#Exercise: if annual revenue is 14000000, and Capital expenditure is 150000000, then please calculate Net present value. Should we invest this project?

npv=calc_NPV(annual_revenue = 8090000000,lifetime_yrs=25, CAPEX=42610000000  )
ifelse(npv>0, "Support","obeject" )

cash_flows_result <- calc_NPV(annual_revenue = 8090000000, lifetime_yrs = 25, CAPEX = 42610000000)
print(cash_flows_result)
#############################
calc_NPV <- function(annual_revenue, i=0.05, lifetime_yrs, CAPEX, OPEX=0){
  costs_op <- rep(OPEX, lifetime_yrs) #operating cost
  revenue <- rep(annual_revenue, lifetime_yrs) 
  t <- seq(1, lifetime_yrs, 1) #output: 1, 2, 3, ...25
  
  NPV <- rep(0, lifetime_yrs)  # Initialize vector to store NPV for each year
  
  for (year in 1:lifetime_yrs) {
    NPV[year] <- (revenue - costs_op)/(1 + i)**year
  }
  
  NPV <- round(NPV, 0)
  return(NPV)
}

npv_vector <- calc_NPV(annual_revenue = 8090000000,lifetime_yrs=25, CAPEX=42610000000)

# Print NPV for each year
for (year in 1:length(npv_vector)) {
  cat("Year", year, ": Flow in =", npv_vector[year], "\n")
}

#2 Levelized cost of electricity (LCOE)=====
#Life_span_generation_kWH is one of the required inputs to estimate the Levelized
#cost of electricity (following function)
Life_span_generation_kWH <- function (yearly_generation_kWH, discount = 0.08, lifetime_yrs = 25){
  t<- seq(1, lifetime_yrs, 1)
  L_S_G <- sum(yearly_generation_kWH/(1+discount)**t)
  return (round(L_S_G,0))
}
generation=Life_span_generation_kWH(yearly_generation_kWH=78550000000)

#NPV of cost. 
calc_LCOE <- function(NPV,Life_span_generation){
  lcoe <- NPV/Life_span_generation
  return(round(lcoe,2))
}


LCOE=calc_LCOE(NPV=71410011539,Life_span_generation=generation)
###############################