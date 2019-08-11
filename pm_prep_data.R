# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script prepares a long-term mark recapture data set for analysis and 
# selects five species for an investigation of partial migration. The only
# requirement is to set the working directory, which should include a 
# subfolder called: 'data'

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Read in raw data
raw.data <- read.table("data/oritag_single_recaps.csv", 
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       stringsAsFactors = FALSE,
                       head = TRUE)

# Create new data frame for exploration, formatting and cleaning
prep.data <- raw.data

# EXPLORE DATA ------------------------------------------------------------

names(prep.data)
head(prep.data)
tail(prep.data)
summary(prep.data)
str(prep.data)

# FORMAT DATA -------------------------------------------------------------

# Coerce tag_date and recapture_date to dates
prep.data$tag_date <- as.POSIXlt(prep.data$tag_date)
prep.data$recapture_date <- as.POSIXlt(prep.data$recapture_date)

# Add tag year column
prep.data$tag_year <- 1900 + as.POSIXlt(prep.data$tag_date)$year
prep.data$tag_year <- as.factor(prep.data$tag_year)  # Convert year to factor

# Add recapture year column
prep.data$recapture_year <- 1900 + as.POSIXlt(prep.data$recapture_date)$year
prep.data$recapture_year <- as.factor(prep.data$recapture_year)  # Convert year to factor

# Add rounded release locality
prep.data$release_locality_r <- round(prep.data$release_locality / 100) * 100

# Direction of migratory behaviour
prep.data$direction <- prep.data$release_locality - prep.data$recapture_locality

# SPECIES SELECTION -------------------------------------------------------

# Create new data frame for species selection
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

# Filter data leaving only top 5
dat <- dat[dat$species == 'LRVS'|
             dat$species == 'RTSH'|
             dat$species == 'SGSH'|
             dat$species == 'GLJN'|
             dat$species == 'SSNP',]

# Drop latent factor levels
dat <- droplevels(dat)

# EXPORT FINAL DATA SET FOR ANALYSIS --------------------------------------

write.table(dat, 
            "data/pm_dat.csv",
            sep = ",",
            row.name = FALSE,
            col.names = TRUE)

# END OF FILE -------------------------------------------------------------
