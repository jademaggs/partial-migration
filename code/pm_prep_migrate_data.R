# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script takes a subset of the 'pm_dat.csv' data set, called
# 'migrate_data' and adds variables for investigating direction of wide-
# ranging movement. The only requirement is to set the working directory, 
# which should include a subfolder called: 'data', with a data set called
# 'migrate_data.csv'.

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Import libraries
library(fossil)  # Calculate bearing between coordinate pairs

# Read in datasets
migrate.dat <- read.table("data/migrate_data.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       head = TRUE)

# ADD DIRECTION VARIABLES TO DATA SET -------------------------------------

# Add polar bearing column to dataframe
migrate.dat$bearing.deg <-
  earth.bear(
    migrate.dat$release_long,
    migrate.dat$release_lat,
    migrate.dat$recap_long,
    migrate.dat$recap_lat
  )

# Convert bearing in degrees to radians
migrate.dat$bearing <- migrate.dat$bearing * (pi/180)

# Define function to convert polar bearing to cartesian coordinates
polar2cart<-function(x,y, dist,bearing,as.deg=FALSE){
  ## Translate Polar coordinates into Cartesian coordinates
  ## based on starting location, distance, and bearing
  ## as.deg indicates if the bearing is in degrees (T) or radians (F)
  
  if(as.deg){
    ##if bearing is in degrees, convert to radians
    bearing=bearing*pi/180
  }
  
  newx<-x+dist*sin(bearing)  ##X
  newy<-y+dist*cos(bearing)  ##Y
  return(list("x"=newx,"y"=newy))
}

# Convert polar bearing to cartesian coordinates
migrate.dat$cartesian_x <- polar2cart(0,0, 
                                    migrate.dat$km_moved,
                                    migrate.dat$bearing.deg, TRUE)$x

migrate.dat$cartesian_y <- polar2cart(0,0, 
                                      migrate.dat$km_moved,
                                      migrate.dat$bearing.deg, TRUE)$y

# Add azimuths
migrate.dat$azimuth_y <- sin(migrate.dat$bearing)
migrate.dat$azimuth_x <- cos(migrate.dat$bearing)

# EXPORT FINAL DATA SET FOR ANALYSIS --------------------------------------

write.table(migrate.dat, 
            "data/migrate_direction_data.csv",
            sep = ",",
            row.name = FALSE,
            col.names = TRUE)

# END OF FILE -------------------------------------------------------------
