# AUTHOR:       JQ Maggs
# CONTACT:      jademaggs@gmail.com
# DATE:         2019-08-04

# OVERVIEW ----------------------------------------------------------------

# This script plots the magnitude and direction of wide-ranging movement 
# for the five study species. This is done in terms of bioregion, life-stage, 
# and species. The only requirement is to set the working directory, which 
# should include a sub-folder called 'data', with a data set called
# 'migrate_data.csv'. 

# SETUP -------------------------------------------------------------------

# Erase all previous variables from memory
rm(list = ls())

# Set working directory
setwd("c:/.../")

# Read in dataset
migrate.dat <- read.table("data/migrate_data.csv", 
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

# PLOT MOVEMENT DIRECTION -------------------------------------------------

dat <- raw.data

# Plot direction RTSH ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Carcharias taurus")), family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Delagoa']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Delagoa']), 
       length = 0.08)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'RTSH' & 
                                 migrate.dat$bioregion == 'Natal']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'RTSH' &
                                 migrate.dat$bioregion == 'Natal']), 
       length = 0.08)

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
       length = 0.08)

rtsh_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'RTSH' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'RTSH' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

text(-56, -180, "Delagoa \nn=1", cex = 1.2)
text(-180, -280, "Natal \nn=40, r=0.64", cex = 1.2)
text(200, 206, "Agulhas \nn=31, r=0.15", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')

# Plot direction LRVS ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lichia amia")), family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'LRVS' & 
                                      migrate.dat$bioregion == 'Natal']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'LRVS' &
                                      migrate.dat$bioregion == 'Natal']), 
       length = 0.08)

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
       length = 0.08)

lrvs_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'LRVS' &
                              migrate.dat$bioregion == 'Agulhas']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'LRVS' &
                                migrate.dat$bioregion == 'Agulhas']) ^ 2)

text(-174, -180, "Natal \nn=23, r=0.01", cex = 1.2)
text(230, 280, "Agulhas\nn=95, r=0.44", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')

# Plot direction SGSH ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Triakis megalopterus")), family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SGSH' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SGSH' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.08)


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
       length = 0.08)

sgsh_namib_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'SGSH' &
                              migrate.dat$bioregion == 'Namib']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'SGSH' &
                                migrate.dat$bioregion == 'Namib']) ^ 2)

text(60, 30, "Agulhas\nn=16, r=0.21", cex = 1.2)
text(-20, 20, "Namib\nn=7, r=0.15", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')

# Plot direction GLJN ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-500, 500), ylim = c(-500, 500), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Dichistius capensis")), family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'GLJN' & 
                                      migrate.dat$bioregion == 'Agulhas']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'GLJN' &
                                      migrate.dat$bioregion == 'Agulhas']), 
       length = 0.08)

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
       length = 0.08)

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
       length = 0.08)

gljn_agulhas_r <-  # Root mean square of x and y azimuths
  sqrt(
    mean(
      migrate.dat$azimuth_y[migrate.dat$species == 'GLJN' &
                              migrate.dat$bioregion == 'Namib']) ^ 2 
    +
      mean(
        migrate.dat$azimuth_x[migrate.dat$species == 'GLJN' &
                                migrate.dat$bioregion == 'Namib']) ^ 2)

text(120, 100, "Agulhas\nn=232, r=0.05", cex = 1.2)
text(250, -100, "Namaqua\nn=46, r=0.87", cex = 1.2)
text(250, -460, "Namib\nn=6, r=0.66", cex = 1.2)
text(0, 480, "North", cex = 1.2, family='D')

# Plot direction SSNP ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lutjanus rivulatus")), family='D')

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SSNP' & 
                                      migrate.dat$bioregion == 'Delagoa']), 
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Delagoa']), 
       length = 0.08)

arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Natal']),
       mean(migrate.dat$cartesian_y[migrate.dat$species == 'SSNP' &
                                      migrate.dat$bioregion == 'Natal']),
       length = 0.08)

text(35, 55, "Delagoa\nn=1", cex = 1.2)
text(-7, 55, "Natal\nn=1", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')

# Plot bioregion ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main="All species")

# Add Delagoa bioregion
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$bioregion == 'Delagoa']),
       mean(migrate.dat$cartesian_y[migrate.dat$bioregion == 'Delagoa']),
       length = 0.08)

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
       length = 0.08)

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
       length = 0.08)

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
       length = 0.08)

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
       length = 0.08)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$bioregion == 'Namib']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$bioregion == 'Namib']) ^ 2)

text(-40, -140, "Delagoa\nn=2, r=0.05", cex = 1.2)
text(-200, -250, "Natal\nn=65, r=0.38", cex = 1.2)
text(200, 150, "Agulhas\nn=374, r=0.14", cex = 1.2)
text(220, -70, "Namaqua\nn=46, r=0.87", cex = 1.2)
text(100, -250, "Namib\nn=13, r=0.38", cex = 1.2)
text(0, 300, "North", cex = 1.2, family = 'D')

# Plot life_stage (LRVS) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-300, 300), ylim = c(-300, 300), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lichia amia")), family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'LRVS_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'LRVS_ad']),
       length = 0.08)

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
       length = 0.08)

# Root mean square of x and y azimuths
sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'LRVS_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'LRVS_juv']) ^ 2)

text(220, 45, "Juveniles\nn=42, r=0.38", cex = 1.2)
text(180, 220, "Adults\nn=76, r=0.34", cex = 1.2)
text(0, 300, "North", cex = 1.2, family='D')

# Plot life_stage (RTSH) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-200, 200), ylim = c(-200, 200), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Carcharias taurus")), family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'RTSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'RTSH_ad']),
       length = 0.08)

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
       length = 0.08)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'RTSH_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'RTSH_juv']) ^ 2)

text(-100, -160, "Adults\nn=57, r=0.36", cex = 1.2)
text(80,-30, "Juveniles\nn=15, r=0.09", cex = 1.2)
text(0, 200, "North", cex = 1.2, family='D')

# Plot life_stage (SGSH) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Triakis megalopterus")), family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SGSH_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SGSH_ad']),
       length = 0.08)

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
       length = 0.08)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'SGSH_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'SGSH_juv']) ^ 2)


text(70, 30, "Adults\nn=18, r=0.11", cex = 1.2)
text(30,-45, "Juveniles\nn=5, r=0.37", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')

# Plot life_stage (GLJN) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Dichistius capensis")), family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'GLJN_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'GLJN_ad']),
       length = 0.08)

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
       length = 0.08)

sqrt(
  mean(
    migrate.dat$azimuth_y[migrate.dat$species_stage == 'GLJN_juv']) ^ 2 
  +
    mean(
      migrate.dat$azimuth_x[migrate.dat$species_stage == 'GLJN_juv']) ^ 2)

text(65, 30, "Adults\nn=272, r=0.17", cex = 1.2)
text(65,-20, "Juveniles\nn=13, r=0.41", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')

# Plot life_stage (SSNP) ----
plot(1, type='n', axes = TRUE, xlab="km", ylab="km",
     xlim = c(-100, 100), ylim = c(-100, 100), cex.axis = 1.2, cex.lab=1.2, las=1)
abline(v=c(0,0), col='grey')
abline(h=c(0,0), col='grey')
title(main=substitute(italic("Lutjanus rivulatus")), family = 'D')

# Add adults
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_ad']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_ad']),
       length = 0.08)

# Add juveniles
arrows(0,0,
       mean(migrate.dat$cartesian_x[migrate.dat$species_stage == 'SSNP_juv']),
       mean(migrate.dat$cartesian_y[migrate.dat$species_stage == 'SSNP_juv']),
       length = 0.08)

text(-10, 60, "Adults\nn=1", cex = 1.2)
text(40,60, "Juveniles\nn=1", cex = 1.2)
text(0, 100, "North", cex = 1.2, family='D')

# END OF FILE -------------------------------------------------------------
