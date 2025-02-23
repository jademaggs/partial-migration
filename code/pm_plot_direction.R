# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script plots the magnitude and direction of wide-ranging movement 
# for the five study species. This is done in terms of species, bioregion 
# and life-stage The only requirement is to set the working directory, which 
# should include a sub-folder called 'data', with a data set called
# 'migrate_data.csv'. All plots are exported to 'output' sub-folder.  

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Read in dataset
migrate.dat <- read.table("data/migrate_direction_data.csv", 
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

# PLOT MOVEMENT DIRECTION ACCORDING TO SPECIES ---------------------------

species <- list('RTSH','LRVS','SGSH','GLJN','SSNP')
bioregion <- list('Namib','Namaqua','Agulhas','Natal','Delagoa')

# Plot direction RTSH ----
tiff(file='output/direct_rtsh.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Carcharias taurus")), family='D')

# Add direction vectors for each bioregion
for (i in bioregion){
  arrows(0,0,
         mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                        migrate.dat$bioregion == i]), 
         mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' & 
                                        migrate.dat$bioregion == i]), 
         length = 0.08) 
}

text(-56, -180, "Delagoa \nn=1", cex = 1.2)
text(-180, -280, "Natal \nn=40, r=0.64", cex = 1.2)
text(200, 206, "Agulhas \nn=31, r=0.15", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')
dev.off()

# Plot direction LRVS ----
tiff(file='output/direct_lrvs.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lichia amia")), family='D')

# Add direction vectors for each bioregion
for (i in bioregion){
  arrows(0,0,
         mean(migrate.dat$cartesian_x[migrate.dat$species == 'LRVS' & 
                                        migrate.dat$bioregion == i]), 
         mean(migrate.dat$cartesian_y[migrate.dat$species == 'LRVS' & 
                                        migrate.dat$bioregion == i]), 
         length = 0.08) 
}

text(-174, -180, "Natal \nn=23, r=0.01", cex = 1.2)
text(230, 280, "Agulhas\nn=95, r=0.44", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')
dev.off()

# Plot direction SGSH ----
tiff(file='output/direct_sgsh.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), 
     cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Triakis megalopterus")), family='D')

# Add direction vector for Agulhas Bioregion
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SGSH' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SGSH' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.08)

# Add direction vector for Namib Bioregion
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SGSH' & 
                                      migrate.dat$bioregion == 'Namib']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SGSH' &
                                      migrate.dat$bioregion == 'Namib']), 
       length = 0.08)

text(60, 30, "Agulhas\nn=16, r=0.21", cex = 1.2)
text(-20, 20, "Namib\nn=7, r=0.15", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')
dev.off()

# Plot direction GLJN ----
tiff(file='output/direct_gljn.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-500, 500), ylim = c(-500, 500), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Dichistius capensis")), family='D')

# Add direction vectors for each bioregion
for (i in bioregion){
  arrows(0,0,
         mean(migrate.dat$cartesian_x[migrate.dat$species == 'GLJN' & 
                                        migrate.dat$bioregion == i]), 
         mean(migrate.dat$cartesian_y[migrate.dat$species == 'GLJN' & 
                                        migrate.dat$bioregion == i]), 
         length = 0.08) 
}

text(120, 100, "Agulhas\nn=232, r=0.05", cex = 1.2)
text(250, -100, "Namaqua\nn=46, r=0.87", cex = 1.2)
text(250, -460, "Namib\nn=6, r=0.66", cex = 1.2)
text(0, 480, "North", cex = 1.2, family='D')
dev.off()

# Plot direction SSNP ----
tiff(file='output/direct_ssnp.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lutjanus rivulatus")), family='D')

# Add direction vectors for each bioregion
for (i in bioregion){
  arrows(0,0,
         mean(migrate.dat$cartesian_x[migrate.dat$species == 'SSNP' & 
                                        migrate.dat$bioregion == i]), 
         mean(migrate.dat$cartesian_y[migrate.dat$species == 'SSNP' & 
                                        migrate.dat$bioregion == i]), 
         length = 0.08) 
}

text(35, 55, "Delagoa\nn=1", cex = 1.2)
text(-7, 55, "Natal\nn=1", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')
dev.off()

# PLOT MOVEMENT DIRECTION ACCORDING TO BIOREGION -------------------------
tiff(file='output/direct_bioreg.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="All species")

# Add direction vectors for each bioregion
for (i in bioregion){
  arrows(0,0,
         mean(migrate.dat$cartesian_x[migrate.dat$bioregion == i]),
         mean(migrate.dat$cartesian_y[migrate.dat$bioregion == i]), 
         length = 0.08) 
}

text(-40, -140, "Delagoa\nn=2, r=0.05", cex = 1.2)
text(-200, -250, "Natal\nn=65, r=0.38", cex = 1.2)
text(200, 150, "Agulhas\nn=374, r=0.14", cex = 1.2)
text(220, -70, "Namaqua\nn=46, r=0.87", cex = 1.2)
text(100, -250, "Namib\nn=13, r=0.38", cex = 1.2)
text(0, 300, "North", cex = 1.2, family = 'D')
dev.off()

# PLOT MOVEMENT DIRECTION ACCORDING TO LIFE-STAGE -------------------------

# Plot life_stage (LRVS) ----
tiff(file='output/direct_lrvs_life_stage.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lichia amia")), family = 'D')

# Add direction vector for adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'LRVS_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'LRVS_ad']),
       length = 0.08)

# Add direction vector for juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'LRVS_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'LRVS_juv']),
       length = 0.08)

text(220, 45, "Juveniles\nn=42, r=0.38", cex = 1.2)
text(180, 220, "Adults\nn=76, r=0.34", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')
dev.off()

# Plot life_stage (RTSH) ----
tiff(file='output/direct_rtsh_life_stage.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-200, 200), ylim = c(-200, 200), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Carcharias taurus")), family = 'D')

# Add direction vector for adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'RTSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'RTSH_ad']),
       length = 0.08)

# Add direction vector for juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'RTSH_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'RTSH_juv']),
       length = 0.08)

text(-100, -160, "Adults\nn=57, r=0.36", cex = 1.2)
text(80,-30, "Juveniles\nn=15, r=0.09", cex = 1.2)
text(0, 200, "North", cex = 1.2, family='D')
dev.off()

# Plot life_stage (SGSH) ----
tiff(file='output/direct_sgsh_life_stage.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Triakis megalopterus")), family = 'D')

# Add direction vector for adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SGSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SGSH_ad']),
       length = 0.08)

# Add direction vector for juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SGSH_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SGSH_juv']),
       length = 0.08)

text(70, 30, "Adults\nn=18, r=0.11", cex = 1.2)
text(30,-45, "Juveniles\nn=5, r=0.37", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')
dev.off()

# Plot life_stage (GLJN) ----
tiff(file='output/direct_gljn_life_stage.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Dichistius capensis")), family = 'D')

# Add direction vector for adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'GLJN_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'GLJN_ad']),
       length = 0.08)

# Add direction vector for juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'GLJN_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'GLJN_juv']),
       length = 0.08)

text(65, 30, "Adults\nn=272, r=0.17", cex = 1.2)
text(65,-20, "Juveniles\nn=13, r=0.41", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')
dev.off()

# Plot life_stage (SSNP) ----
tiff(file='output/direct_ssnp_life_stage.tiff')
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), 
     cex.axis = 1.2, cex.lab=1.2, las=1)

abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lutjanus rivulatus")), family = 'D')

# Add direction vector for adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_ad']),
       length = 0.08)

# Add direction vector for juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_juv']),
       length = 0.08)

text(-10, 60, "Adults\nn=1", cex = 1.2)
text(40,60, "Juveniles\nn=1", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')
dev.off()

# END OF FILE -------------------------------------------------------------
