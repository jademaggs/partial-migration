# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script first develops a binomial model of fish movement using a 
# long-term mark-recapture data set to investigate the probability of 
# wide-ranging movement. The script then develops a Gaussian model of 
# wide-ranging movement to investigate direction of movement. The only
# requirement is to set the working directory, which should include a 
# sub-folder called 'data', containing a data set called 'pm_dat.csv'.

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Import libraries
library(lmtest)  # Likelihood ratio test
library(fossil)  # Calculate bearing between coordinate pairs
library(vegan)  # MDS analyses
library(multcomp)  # Posthoc testing for GLM

# Read in datasets
raw.data <- read.table("data/pm_dat.csv", # From pm_prep_data.R
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       stringsAsFactors = FALSE,
                       head = TRUE)

migrate.dat <- read.table("data/migrate_data.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       head = TRUE)

# DEVELOP BINOMIAL MODEL OF MOVEMENT --------------------------------------

dat <- raw.data

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

# Predict bioregion
prediction_bioregion <- 
  data.frame(species = rep(c("GLJN","GLJN","GLJN","GLJN","GLJN"), each=1),
             bioregion = rep(c("Delagoa","Natal","Agulhas","Namaqua","Namib"),times=1),
             life_stage = rep("AD", each=5))
prediction_bioregion <- cbind(prediction_bioregion, predict(m3, prediction_bioregion, type="response"))

# DEVELOP GAUSSIAN MODEL OF MOVEMENT DIRECTION ----------------------------

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

# END OF FILE -------------------------------------------------------------
