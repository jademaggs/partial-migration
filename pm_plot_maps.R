# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script plots maps showing the spatial distribution of resident and 
# wide-ranging behaviour for each of the five study species and two 
# life-stages. The only requirement is to set the working directory, which 
# should include a sub-folder called 'data', with a data set called
# 'map_data.csv'. All plots are exported to 'output' sub-folder.  

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Import libraries
library(maps)  # ver3.3.0 - Plotting maps
library(mapdata)  # Map data sources 

# Read in datasets
map.data <- read.table("data/map_data.csv",
                       sep=",",
                       dec=".",
                       na.strings="NA",
                       strip.white=TRUE,
                       head = TRUE)

# Import windows fonts
windowsFonts(A=windowsFont("Arial Black"),
             B=windowsFont("Arial Narrow"),
             C=windowsFont("Arial Regular"),
             D=windowsFont("Arial Italic"))

# PLOT RESIDENT VS MIGRATORY BEHAVIOUR SPATIALLY --------------------------

long <- map.data$longitude
lat <- map.data$latitude

# RTSH - Carcharias taurus ----

# Juveniles
tiff(file='output/behaviour_rtsh_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'RTSH_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'RTSH_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'RTSH_juv' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'RTSH_juv') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text(substitute(italic("Carcharias taurus")), x=23, y=-25, family="D")
text("Juveniles", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# Adults
tiff(file='output/behaviour_rtsh_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'RTSH_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'RTSH_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'RTSH_ad' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'RTSH_ad') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text("Adults", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# LRVS - Lichia amia ----

# Juvenile
tiff(file='output/behaviour_lrvs_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'LRVS_juv' & 
              map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'LRVS_juv'& 
             map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'LRVS_juv' & 
              map.data$move_type == 'resident'],
     lat[(map.data$species_stage == 'LRVS_juv') & 
           map.data$move_type == 'resident'], 
     pch="*", col='black', cex=5)

text(substitute(italic("Lichia amia")), x=23, y=-25, family="D")
text("Juveniles", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# Adults
tiff(file='output/behaviour_lrvs_ad.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'LRVS_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'LRVS_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'LRVS_ad' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'LRVS_ad') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text("Adults", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# SGSH - Triakis meglaopterus ----

# Juveniles
tiff(file='output/behaviour_sgsh_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'SGSH_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SGSH_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'SGSH_juv' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'SGSH_juv') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text(substitute(italic("Triakis megalopterus")), x=23, y=-25, family="D")
text("Juveniles", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# Adults
tiff(file='output/behaviour_sgsh_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'SGSH_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SGSH_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'SGSH_ad' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'SGSH_ad') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text("Adults", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# GLJN - Dichistius capensis ----

# Juveniles
tiff(file='output/behaviour_gljn_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'GLJN_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'GLJN_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'GLJN_juv' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'GLJN_juv') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text(substitute(italic("Dichistius capensis")), x=23, y=-25, family="D")
text("Juveniles", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# Adults
tiff(file='output/behaviour_gljn_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'GLJN_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'GLJN_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'GLJN_ad' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'GLJN_ad') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text("Adults", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# SSNP - Lutjanus rivulatus ----

# Juveniles
tiff(file='output/behaviour_ssnp_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'SSNP_juv' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SSNP_juv'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'SSNP_juv' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'SSNP_juv') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text(substitute(italic("Lutjanus rivulatus")), x=23, y=-25, family="D")
text("Juveniles", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# Adults
tiff(file='output/behaviour_ssnp_juv.tiff')
map("worldHires", interior=FALSE, fill=TRUE, col='white', border='white', 
    xlim=c(11,38), ylim=c(-36, -18), mar=c(1,1,1,1))

points(long[map.data$species_stage == 'SSNP_ad' & map.data$move_type == 'migratory'],
       lat[map.data$species_stage == 'SSNP_ad'& map.data$move_type == 'migratory'], 
       pch=2,col='grey39', cex=4, lwd =5)

points(long[map.data$species_stage == 'SSNP_ad' & map.data$move_type == 'resident'],
       lat[(map.data$species_stage == 'SSNP_ad') & map.data$move_type == 'resident'], 
       pch="*", col='black', cex=5)

text("Adults", x=23, y=-32, family="C")
map('world', interior = F, add = T)
dev.off()

# END OF FILE -------------------------------------------------------------
