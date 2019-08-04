# AUTHOR:       Jade Q. Maggs MSc
# POSITION:     Assistant Scientist
# AFFILIATION:  Oceanographic Research Institute / University of KwaZulu-Natal
# CONTACT:      email: jmaggs@ori.org.za / mobile: +2783 515 1079 

# 1. SCRIPT METADATA ------------------------------------------------------

# INSTRUCTIONS    Execute "SETUP" first then individual sections in "MAIN".
# OVERVIEW:       Analyse ORITag dat for evidence of partial migration (PhD)
# START DATE:     2016-04-11
# INPUTS:         Imports text file from "/data" folder in working directory
# OUTPUTS:        Exports all output to "/output" folder in working directory
# DATA ORIGIN:    ORITag
# DATE EXTRACTED: 2015-05-04 12:04
# DATA COVERAGE:  Southern Africa
#                 1984-01-01 - 2015-05-04
#                 marine / estuarine
#                 multiple recaptures excluded
#                 --

# 2. SETUP ----------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("C:/Users/jmaggs/Dropbox/PhD/thesis/ch5_partial_migration/")

# 2.1. IMPORT OTHER SOURCE CODE

# 2.2. IMPORT PACKAGES
library(ggplot2)
library(scales)
# library(spatial)
# library(MASS)  # negative binomial glm
library(maps)  # Plotting maps
library(mapdata)  # Map data sources 
library(ggmap)  # ggplot maps
library(lmtest)  # Likelihood ratio test
library(tweedie)
library(fossil)  # Calculate bearing between coordinate pairs
library(vegan)
library(multcomp)  # Posthoc testing for GLM


# 2.3. FUNCTION DEFINITIONS

# Convert polar bearing to cartesian coordinates
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

# 2.4. IMPORT DATASETS

# Read in dataset
raw.data <- read.table("data/oritag_single_recaps.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       stringsAsFactors = FALSE,
                       head = TRUE)

map.data <- read.table("data/map_data.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       head = TRUE)

migrate.dat <- read.table("data/migrate_data.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       head = TRUE)

# Create new dataset for exploration, formatting and cleaning
prep.data <- raw.data

# 2.5. EXPLORE AND CLEAN data
names(prep.data)
head(prep.data)
tail(prep.data)
summary(prep.data)
str(prep.data)

# Species list according to no.of recaptures
aggregate(prep.data$tag_number,list(prep.data$species),length)

# Look for zero days free with km greater than 0
prep.data[prep.data$days_free == 0 & prep.data$km_moved > 0, ]

# 2.6. FORMAT data

# Coerce tag_date and recapture_date to dates
prep.data$tag_date <- as.POSIXlt(prep.data$tag_date)
prep.data$recapture_date <- as.POSIXlt(prep.data$recapture_date)

# Add tag year column
prep.data$tag_year <- 1900 + as.POSIXlt(prep.data$tag_date)$year
prep.data$tag_year <- as.factor(prep.data$tag_year)  # Convert year to factor

# Add recapture year column
prep.data$recapture_year <- 1900 + as.POSIXlt(prep.data$recapture_date)$year
prep.data$recapture_year <- as.factor(prep.data$recapture_year)  # Convert year to factor

# Remove unecessary data
prep.data <- 
  data.frame(prep.data[prep.data$tag_date >= 1984-01-01, ])# Retain only 1984-2015

prep.data <- 
  data.frame(prep.data[prep.data$recapture_locality != 9999, ])# Drop unknown locality

# 2.7. GLOBAL DECLARATIONS

# Import windows fonts
windowsFonts(A=windowsFont("Arial Black"),
             B=windowsFont("Arial Narrow"),
             C=windowsFont("Arial Regular"),
             D=windowsFont("Arial Italic"))

# 3. MAIN -----------------------------------------------------------------

# Initial spp selection (Table 1) ----
dat <- prep.data

# Create new binary variable 'move_type' - resident/migratory
dat$move_type <- 'unknown'
dat$move_type[dat$km_moved <= 5 & dat$days_free > 365] = 'resident'
dat$move_type[dat$km_moved >= 51 & dat$days_free <= 365] = 'migratory'

# Select only species with known movement types
dat <- dat[dat$move_type != 'unknown',]
length(dat$tag_number)
tapply(dat$species, list(dat$species_stage, dat$move_type), length)
tapply(dat$cluster, list(dat$cluster, dat$move_type), length)

# Resident observations by life stage
write.table(
  tapply(dat$species[dat$move_type == 'resident'],
         list(dat$species[dat$move_type == 'resident'],
              dat$life_stage[dat$move_type == 'resident']),
         length),
  file = 'clipboard',
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE
)

# Migratory observations by life stage
write.table(
  tapply(dat$species[dat$move_type == 'migratory'], 
       list(dat$species[dat$move_type == 'migratory'], 
            dat$life_stage[dat$move_type == 'migratory']), 
       length),
  file = 'clipboard',
  sep = '\t',
  row.names = TRUE,
  col.names = TRUE
)

# Add some variables ------------------------------------------------------

# Add rounded release locality
dat$release_locality_r <- round(dat$release_locality / 100) * 100

# Direction of migratory behaviour
dat$direction <- dat$release_locality - dat$recapture_locality

# Filter data leaving only top 5 ----

dat <- dat[dat$species == 'LRVS'|
             dat$species == 'RTSH'|
             dat$species == 'SGSH'|
             dat$species == 'GLJN'|
             dat$species == 'SSNP',]

dat <- droplevels(dat)

# Plot resident vs migratory behaviour spatially using XY axes ----
ggplot(dat[dat$cluster == 'I',], 
       aes(release_locality_r, species_scientific_stage)) + 
  xlab('Tag-release locality') +
  ylab('Species') +
  xlim(8000,2000) + 
  theme_bw() +
  theme(text = element_text(size=10, family="C")) +
  geom_point(data=dat[dat$move_type == 'resident' & dat$cluster == 'I',],
             col='red',shape=16,size=1) +
  geom_point(data=dat[dat$move_type == 'migratory' & dat$cluster == 'I',],
             col='darkgreen',shape=2,size=2) +
  annotate("text", x=2500, y='Triakis megalopterus_AD', 
           label="Type I\nWide ranging\nspecies",
           size=3, family="C")

# Plot resident vs migratory behaviour spatially using maps ----
long <- map.data$longitude
lat <- map.data$latitude

# LRVS juv
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'LRVS_juv' & map.data$move_type == 'resident'],
     lat[map.data$species_stage == 'LRVS_juv'& map.data$move_type == 'resident'], 
     xlim=c(13,38), ylim=c(-36, -22),
     xlab="Longitude",
     ylab="Latitutde",
     pch=20, col='red')
points(long[map.data$species_stage == 'LRVS_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'LRVS_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Lichia amia", x=33, y=-33, family="D")
text("(Juvenile)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# LRVS ad
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'LRVS_ad' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'LRVS_ad'& map.data$move_type == 'resident'], 
       xlim=c(13,38), ylim=c(-36, -22),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'LRVS_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'LRVS_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Lichia amia", x=33, y=-33, family="D")
text("(Adult)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# RTSH juv
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'RTSH_juv' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'RTSH_juv'& map.data$move_type == 'resident'], 
       xlim=c(13,38), ylim=c(-36, -22),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'RTSH_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'RTSH_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Carcharias taurus", x=33, y=-33, family="D")
text("(Juvenile)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# RTSH ad
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'RTSH_ad' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'RTSH_ad'& map.data$move_type == 'resident'], 
       xlim=c(13,38), ylim=c(-36, -22),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'RTSH_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'RTSH_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Carcharias taurus", x=33, y=-33, family="D")
text("(Adult)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# SSNP juv
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'SSNP_juv' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'SSNP_juv'& map.data$move_type == 'resident'], 
       xlim=c(13,38), ylim=c(-36, -22),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'SSNP_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SSNP_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Lutjanus rivulatus", x=33, y=-33, family="D")
text("(Juvenile)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# SSNP ad
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(13,38), ylim=c(-36, -22), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'SSNP_ad' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'SSNP_ad'& map.data$move_type == 'resident'], 
       xlim=c(13,38), ylim=c(-36, -22),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'SSNP_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SSNP_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Lutjanus rivulatus", x=33, y=-33, family="D")
text("(Adult)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# GLJN juv
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(12,38), ylim=c(-36, -18), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'GLJN_juv' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'GLJN_juv'& map.data$move_type == 'resident'], 
       xlim=c(12,38), ylim=c(-36, -18),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'GLJN_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'GLJN_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Dichistius capensis", x=33, y=-33, family="D")
text("(Juvenile)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# GLJN ad
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(12,38), ylim=c(-36, -18), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'GLJN_ad' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'GLJN_ad'& map.data$move_type == 'resident'], 
       xlim=c(12,38), ylim=c(-36, -18),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'GLJN_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'GLJN_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Dichistius capensis", x=33, y=-33, family="D")
text("(Adult)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# SGSH juv
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(12,38), ylim=c(-36, -18), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'SGSH_juv' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'SGSH_juv'& map.data$move_type == 'resident'], 
       xlim=c(12,38), ylim=c(-36, -18),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'SGSH_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SGSH_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Triakis megalopterus", x=33, y=-33, family="D")
text("(Juvenile)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# SGSH ad
map("worldHires", interior=FALSE, fill=TRUE, col="white",
    xlim=c(12,38), ylim=c(-36, -18), mar=c(1,1,1,1))
points(long[map.data$species_stage == 'SGSH_ad' & map.data$move_type == 'resident'],
       lat[map.data$species_stage == 'SGSH_ad'& map.data$move_type == 'resident'], 
       xlim=c(12,38), ylim=c(-36, -18),
       xlab="Longitude",
       ylab="Latitutde",
       pch=20, col='red')
points(long[map.data$species_stage == 'GLJN_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'GLJN_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='green', cex=2)
text("Triakis megalopterus", x=33, y=-33, family="D")
text("(Adult)", x=33, y=-34, family="C")
text("SOUTH \nAFRICA", x=23, y=-30, family="C")

# Plot resident vs migratory behaviour temporally (Type I) ----
ggplot(dat[dat$cluster == 'I',], 
       aes(tag_year, species_scientific_stage)) + 
  xlab('Year') +
  ylab('Species') +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(text = element_text(size=10, family="C")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_point(data=dat[dat$move_type == 'resident' & dat$cluster == 'I',],
             col='red',shape=16,size=1) +
  geom_point(data=dat[dat$move_type == 'migratory' & dat$cluster == 'I',],
             col='green',shape=2,size=2) +
  annotate("text", x=28, y='Triakis megalopterus_AD', 
           label="Type I\nWide-ranging species",
           size=3, family="C") 

# Plot resident vs migratory behaviour temporally (Type IIa) ----
ggplot(dat, 
       aes(tag_year, species_scientific_stage)) + 
  geom_point(data=dat[dat$move_type == 'resident' & dat$cluster == 'IIa',],
             aes(tag_year, species_scientific_stage),
             col='red',shape=16,size=1) +
  geom_point(data=dat[dat$move_type == 'migratory' & dat$cluster == 'IIa',],
             aes(tag_year, species_scientific_stage),
             col='green',shape=2,size=2) +
  xlab('Year') +
  ylab('Species') +
  scale_x_discrete(drop=FALSE) +
  theme_bw() +
  theme(text = element_text(size=10, family="C")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  annotate("text", x=6, y='Lutjanus rivulatus_JUV', 
           label="Type II\nResident species",
           size=3, family="C") 

# Plot direction according to distance (Type I) ----
ggplot(dat) + 
  geom_point(data=dat[dat$move_type == 'resident'& dat$cluster == 'I', ],
             aes(direction, species_scientific_stage),
             col='red',shape=16, size=1) +
  geom_point(data=dat[dat$move_type == 'migratory' & dat$cluster == 'I', ],
             aes(direction, species_scientific_stage),
             col='green',shape=2, size=2) +
  xlab('<~~ West  |  Distance moved (km)  |  East ~~>') +
  ylab('Species') +
  xlim(-3000, 3000) + 
  theme_bw() +
  theme(text = element_text(size=10, family="C")) +
  annotate("text", x=2250, y='Triakis megalopterus_AD', 
           label="Type I\nWide-ranging species",
           size=3, family="C") 

# Plot direction according to distance (Type IIa) ----
ggplot(dat) + 
  geom_point(data=dat[dat$move_type == 'resident'& dat$cluster == 'IIa', ],
             aes(direction, species_scientific_stage),
             col='red',shape=16, size=1) +
  geom_point(data=dat[dat$move_type == 'migratory' & dat$cluster == 'IIa', ],
             aes(direction, species_scientific_stage),
             col='green',shape=2, size=2) +
  xlab('<~~ West  |  Distance moved (km)  |  East ~~>') +
  ylab('Species') +
  xlim(-3000, 3000) + 
  theme_bw() +
  theme(text = element_text(size=10, family="C")) +
  annotate("text", x=2250, y='Lutjanus rivulatus_JUV',
           label="Type II\nResident species",
           size=3, family="C") 

# Binomial model of movement behaviour ----
dat$move_type_binomial <- dat$move_type
dat$move_type_binomial[dat$move_type == "resident"] = 0
dat$move_type_binomial[dat$move_type == "migratory"] = 1
dat$move_type_binomial <- as.numeric(dat$move_type_binomial)
dat$species <- as.factor(dat$species)
dat$life_stage <- as.factor(dat$life_stage)
dat$bioregion <- as.factor(dat$bioregion)

test_glm <- glm(move_type_binomial ~ species + 
            life_stage + 
            bioregion,
          family=binomial,
          data=dat)
 
summary(test_glm)
anova(test_glm, test="Chisq")

# Null only
m0 <- glm(move_type_binomial ~ 1,data=dat, family=binomial)
summary(m0)
AIC(m0)
BIC(m0)

# Species
m1 <- glm(move_type_binomial ~ species,data=dat, family=binomial)
summary(m1)
AIC(m1)
BIC(m1)
lrtest(m0, m1) # Does species improve model

# Species + life_stage
m2 <- glm(move_type_binomial ~ species +
            life_stage,data=dat, family=binomial)
summary(m2)
AIC(m2)
BIC(m2)
lrtest(m1, m2) # Does life stage improve model

# Species + life_stage + bioregion
m3 <- glm(move_type_binomial ~ species +
            life_stage + bioregion,data=dat, family=binomial)
summary(m3)
AIC(m3)
BIC(m3)
lrtest(m2, m3) # Does bioregion improve model

# Post hoc testing
summary(glht(m3, mcp(species="Tukey")))
summary(glht(m3, mcp(bioregion="Tukey")))

# Final model summary
write.table(anova(m0, m1, m2, m3), file='clipboard', sep='\t')
anova(m3, test="Chisq")

# Predict juveniles
prediction_juv <- 
  data.frame(species = rep(c("GLJN","LRVS","RTSH","SGSH","SSNP"), each=5),
             bioregion = rep(c("Delagoa","Natal","Agulhas","Namaqua","Namib"),times=5),
             life_stage = rep("JUV", each=25))
             


prediction_juv <- cbind(prediction_juv, predict(m3, prediction_juv, type="response"))
write.table(prediction_juv, file="clipboard", sep="\t",row.names = FALSE)

# predict adults
prediction_ad <- 
  data.frame(species = rep(c("GLJN","LRVS","RTSH","SGSH","SSNP"), each=5),
             bioregion = rep(c("Delagoa","Natal","Agulhas","Namaqua","Namib"),times=5),
             life_stage = rep("AD", each=25))

prediction_ad <- cbind(prediction_ad, predict(m3, prediction_ad, type="response"))
write.table(prediction_ad, file="clipboard", sep="\t",row.names = FALSE)

# Model direction  ----

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

# Run MetaNMDS and extract coefficients for GLM analysis
migrate.matrix <- 
  subset(migrate.dat, select = c('cartesian_x', 'cartesian_y'))
mds1 <- metaMDS(migrate.matrix, distance='euclidean', k=1,trymax=50)
plot(mds1, type='t')
migrate.dat$mds <- mds1$points

# Generalised linear model with nesting
mod <- glm(mds ~ species + bioregion + 
            species/life_stage + bioregion:species/life_stage, 
          family=gaussian, data=migrate.dat)

anova(mod, test='Chisq')

m0 <- glm(mds ~ 1, family=gaussian, data=migrate.dat)
m1 <- glm(mds ~ species, family=gaussian, data=migrate.dat)
m2 <- glm(mds ~ species + bioregion, family=gaussian, data=migrate.dat)
m3 <- glm(mds ~ species + bioregion + species/life_stage, 
          family=gaussian, data=migrate.dat)
m4 <- glm(mds ~ species + bioregion + species/life_stage + 
            bioregion:species/life_stage, family=gaussian, data=migrate.dat)

anova(m0,m1,m2,m3,m4, test='Chisq')

lrtest(m0, m1) # Does species improve model

# Post hoc testing
summary(glht(m2, mcp(species="Tukey")))
summary(glht(m3, mcp(bioregion="Tukey")))

# Plot direction RTSH ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Carcharias taurus", family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Delagoa']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Delagoa']), 
       length = 0.05)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Natal']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' &
                                 migrate.dat$bioregion == 'Natal']), 
       length = 0.05)

rtsh_natal_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'RTSH' &
                                     migrate.dat$bioregion == 'Natal']) ^ 2 
    +
         mean(
           migrate.dat$azimuth_x[migrate.dat$species == 'RTSH' &
                                       migrate.dat$bioregion == 'Natal']) ^ 2)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.05)

rtsh_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'RTSH' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'RTSH' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

text(-56, -180, "Delagoa \nn=1", cex = 0.8)
text(-200, -280, "Natal \nn=40, r=0.64", cex = 0.8)
text(200, 206, "Agulhas \nn=31, r=0.15", cex = 0.8)

text(0, 300, "North", cex = 1, family='D')

# Plot direction LRVS ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Lichia amia", family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'LRVS' & 
                                      migrate.dat$bioregion == 'Natal']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'LRVS' &
                                      migrate.dat$bioregion == 'Natal']), 
       length = 0.05)

lrvs_natal_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'LRVS' &
                              migrate.dat$bioregion == 'Natal']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'LRVS' &
                                migrate.dat$bioregion == 'Natal']) ^ 2)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'LRVS' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'LRVS' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.05)

lrvs_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'LRVS' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'LRVS' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

text(-174, -180, "Natal \nn=23, r=0.01", cex = 0.8)
text(260, 287, "Agulhas\nn=95, r=0.44", cex = 0.8)
text(0, 300, "North", cex = 1, family='D')

# Plot direction SGSH ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Triakis megalopterus", family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SGSH' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SGSH' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.05)


sgsh_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'SGSH' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'SGSH' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SGSH' & 
                                      migrate.dat$bioregion == 'Namib']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SGSH' &
                                      migrate.dat$bioregion == 'Namib']), 
       length = 0.05)

sgsh_namib_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'SGSH' &
                              migrate.dat$bioregion == 'Namib']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'SGSH' &
                                migrate.dat$bioregion == 'Namib']) ^ 2)

text(60, 30, "Agulhas\nn=16, r=0.21", cex = 0.8)
text(-20, 20, "Namib\nn=7, r=0.15", cex = 0.8)
text(0, 100, "North", cex = 1, family='D')

# Plot direction GLJN ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-500, 500), ylim = c(-500, 500))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Dichistius capensis", family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'GLJN' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'GLJN' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.05)

gljn_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'GLJN' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'GLJN' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'GLJN' &
                                      migrate.dat$bioregion == 'Namaqua']),
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'GLJN' &
                                      migrate.dat$bioregion == 'Namaqua']),
       length = 0.05)

gljn_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'GLJN' &
                              migrate.dat$bioregion == 'Namaqua']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'GLJN' &
                                migrate.dat$bioregion == 'Namaqua']) ^ 2)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'GLJN' & 
                                      migrate.dat$bioregion == 'Namib']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'GLJN' &
                                      migrate.dat$bioregion == 'Namib']), 
       length = 0.05)

gljn_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'GLJN' &
                              migrate.dat$bioregion == 'Namib']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'GLJN' &
                                migrate.dat$bioregion == 'Namib']) ^ 2)

text(120, 100, "Agulhas\nn=232 ,r=0.05", cex = 0.8)
text(250, -100, "Namaqua\nn=46, r=0.87", cex = 0.8)
text(250, -460, "Namib\nn=6, r=0.66", cex = 0.8)
text(0, 480, "North", cex = 1, family='D')

# Plot direction SSNP ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Lutjanus rivulatus", family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SSNP' & 
                                      migrate.dat$bioregion == 'Delagoa']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Delagoa']), 
       length = 0.05)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Natal']),
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Natal']),
       length = 0.05)

text(35, 55, "Delagoa\nn=1", cex = 0.8)
text(-7, 55, "Natal\nn=1", cex = 0.8)
text(0, 100, "North", cex = 1, family='D')

# Plot bioregion ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="All species")

# Add Delagoa bioregion
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Delagoa']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Delagoa']),
       length = 0.05)

# Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$bioregion == 'Delagoa']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$bioregion == 'Delagoa']) ^ 2)

# Add Natal
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Natal']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Natal']),
       length = 0.05)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$bioregion == 'Natal']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$bioregion == 'Natal']) ^ 2)

# Add Agulhas
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Agulhas']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Agulhas']),
       length = 0.05)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$bioregion == 'Agulhas']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$bioregion == 'Agulhas']) ^ 2)

# Add Namaqua
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Namaqua']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Namaqua']),
       length = 0.05)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$bioregion == 'Namaqua']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$bioregion == 'Namaqua']) ^ 2)

# Add Namib
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Namib']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Namib']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$bioregion == 'Namib']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$bioregion == 'Namib']) ^ 2)

text(-40, -100, "Delagoa\nn=2, r=0.05", cex = 0.8)
text(-250, -250, "Natal\nn=65, r=0.38", cex = 0.8)
text(200, 100, "Agulhas\nn=374, r=0.14", cex = 0.8)
text(250, -70, "Namaqua\nn=46, r=0.87", cex = 0.8)
text(100, -250, "Namib\nn=13, r=0.38", cex = 0.8)
text(0, 300, "North", cex = 1, family = 'D')

# Plot life_stage (LRVS) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Lichia amia", family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'LRVS_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'LRVS_ad']),
       length = 0.05)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'LRVS_ad']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'LRVS_ad']) ^ 2)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'LRVS_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'LRVS_juv']),
       length = 0.05)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'LRVS_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'LRVS_juv']) ^ 2)

text(250,80, "Juveniles\nn=42, r=0.38", cex = 0.8)
text(180, 210, "Adults\nn=76, r=0.34", cex = 0.8)
text(0, 300, "North", cex = 1, family='D')

# Plot life_stage (RTSH) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-200, 200), ylim = c(-200, 200))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Carcharias taurus", family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'RTSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'RTSH_ad']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'RTSH_ad']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'RTSH_ad']) ^ 2)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'RTSH_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'RTSH_juv']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'RTSH_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'RTSH_juv']) ^ 2)

text(-100, -150, "Adults\nn=57, r=0.36", cex = 0.8)
text(50,-30, "Juveniles\nn=15, r=0.09", cex = 0.8)
text(0, 200, "North", cex = 1, family='D')

# Plot life_stage (SGSH) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Triakis megalopterus", family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SGSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SGSH_ad']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'SGSH_ad']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'SGSH_ad']) ^ 2)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SGSH_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SGSH_juv']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'SGSH_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'SGSH_juv']) ^ 2)


text(70, 15, "Adults\nn=18, r=0.11", cex = 0.8)
text(30,-30, "Juveniles\nn=5, r=0.37", cex = 0.8)
text(0, 100, "North", cex = 1, family='D')

# Plot life_stage (GLJN) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Dichistius capensis", family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'GLJN_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'GLJN_ad']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'GLJN_ad']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'GLJN_ad']) ^ 2)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'GLJN_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'GLJN_juv']),
       length = 0.05)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'GLJN_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'GLJN_juv']) ^ 2)

text(80, 30, "Adults\nn=272, r=0.17", cex = 0.8)
text(80,-20, "Juveniles\nn=13, r=0.41", cex = 0.8)
text(0, 100, "North", cex = 1, family='D')

# Plot life_stage (SSNP) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100))
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="Lutjanus rivulatus", family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_ad']),
       length = 0.05)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_juv']),
       length = 0.05)

text(-10, 60, "Adults\nn=1", cex = 0.8)
text(40,60, "Juveniles\nn=1", cex = 0.8)
text(0, 100, "North", cex = 1, family='D')

# Growth ----
growth.dat <- dat[dat$release_measure == dat$recap_measure &
                    dat$release_length > 0 &
                    dat$recap_length > 0 &
                    dat$recap_length >= dat$release_length &
                    dat$release_measure != 'Unknown',]

growth.dat$growth <- 
  (growth.dat$recap_length - growth.dat$release_length) / growth.dat$days_free

# Mean growth
write.table(
  tapply(
    growth.dat$growth,
    list(growth.dat$species_scientific_stage, growth.dat$move_type),
    median
  ),
  file = 'clipboard',
  sep = '\t',
  row.names = TRUE
)

# Mean release_length
write.table(
  tapply(
    growth.dat$release_length,
    list(growth.dat$species_scientific_stage, growth.dat$move_type),
    median
  ),
  file = 'clipboard',
  sep = '\t',
  row.names = TRUE
)

# Mean recap_length
write.table(
  tapply(
    growth.dat$recap_length,
    list(growth.dat$species_scientific_stage, growth.dat$move_type),
    median
  ),
  file = 'clipboard',
  sep = '\t',
  row.names = TRUE
)

# Number of observations for each species/life-stage
tapply(growth.dat$tag_number, 
       list(growth.dat$species_scientific_stage, growth.dat$move_type),
       length)

# Test of normality
shapiro.test(growth.dat$growth[growth.dat$species_stage == 'RTSH_ad' &
                                 growth.dat$move_type == 'resident'])  # normal

shapiro.test(growth.dat$growth[growth.dat$species_stage == 'RTSH_ad' &
                                 growth.dat$move_type == 'migratory']) 


shapiro.test(growth.dat$growth[growth.dat$species_stage == 'RTSH_juv' &
                                 growth.dat$move_type == 'resident'])  

shapiro.test(growth.dat$growth[growth.dat$species_stage == 'RTSH_juv' &
                                 growth.dat$move_type == 'migratory'])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'RTSH_ad',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'RTSH_juv',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'GLJN_ad',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'GLJN_juv',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'LRVS_ad',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'LRVS_juv',])

wilcox.test(growth ~ move_type, 
            data=growth.dat[growth.dat$species_stage == 'SGSH_ad',])



# 4.END OF FILE -----------------------------------------------------------
