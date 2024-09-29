# Clear R memory
rm(list=ls())
# Close previous graphics windows
graphics.off()
# ---

################################################################################

# LOAD ALL PACKAGES
library("vegan")
library("ggplot2")
library("pairwiseAdonis")
library("ggord")
library("betareg")
library(MuMIn)
# ---

# IMPORT DATA
SALINE <- read.csv("Saline.lagoon.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
SALINE <- na.omit(SALINE)
# ---

# CALCULATE TOTAL ABUNDANCE FOR FISH AND INVERTEBRATES AND ADD TO DATAFRAME
TOTAL.sum <- apply(SALINE[,c("G.spp","S.spp","M.spp","TLGM", "So", "GGM", "FF", "EF", "TSS", "EE","Am","CFL", "EGC", "EMS", "C.spp", "BW", "P")], 1, sum)
TOTAL.sum
SALINE$TOTAL <- TOTAL.sum
FISH.sum <- apply(SALINE[,c("G.spp","M.spp","TLGM", "So", "GGM", "FF", "EF", "TSS", "EE")], 1, sum)
FISH.sum
SALINE$FISH <- FISH.sum
INV.sum <- apply(SALINE[,c("Am","CFL", "EGC", "EMS", "C.spp", "BW", "P")], 1, sum)
INV.sum
SALINE$INV <- INV.sum
# ---

# TRANSFORM TOTAL ABUNDANCES AND ADD TO DATAFRAME (SQUARE ROOT)
T.TOTAL <- sqrt(SALINE$TOTAL)
SALINE$T.TOTAL <- T.TOTAL
T.FISH <- sqrt(SALINE$FISH)
SALINE$T.FISH <- T.FISH
T.INV <- sqrt(SALINE$INV)
SALINE$T.INV <- T.INV
T.SHRIMP <-sqrt(SALINE$S.spp)
SALINE$T.SHRIMP <-T.SHRIMP
# ---

# CALCULATE SHANNON DIVERSITY FOR TOTAL FISH AND INVERTEBRATES AND ADD TO DATAFRAME
total <- SALINE[c("G.spp","S.spp","M.spp","TLGM", "So", "GGM", "FF", "EF", "TSS", "EE","Am","CFL", "EGC", "EMS", "C.spp", "BW", "P")]
total.shan <- diversity(total, index = "shannon")
SALINE$TOTAL.shan <- total.shan  
fish <- SALINE[c("G.spp","M.spp","TLGM", "So", "GGM", "FF", "EF", "TSS", "EE")]
fish.shan <- diversity(fish, index = "shannon")
SALINE$FISH.shan <- fish.shan
INV.data <- SALINE[c("Am","CFL", "EGC", "EMS", "C.spp", "BW", "P")]
INV.shan <- diversity(INV.data, index = "shannon")
SALINE$INV.shan <- INV.shan
# ---

# TRANSFORM AND REPLACE SPECIES ABUNDANCE DATA (HELLINGER)
T.ABUN <- decostand(SALINE[,13:29], method = "hellinger")
SALINE[, 13:29] <- T.ABUN
# ---

# SUBSET DATA
# Sample data subset
SAMPLE <- SALINE[,1:6]
# Standardised water data subset matrix
WATER <- SALINE[,7:12]
WATER.STZ <- scale(WATER)
# Helinger transform water data for PCA
T.WATER <- decostand(WATER, method = "hellinger")
# Transformed fish abundance data subset matrix
T.ABUN <- SALINE[,13:29]
# ---


################################################################################

# ENVIRONMENTAL VARIABILITY

# Create function that tells you min max (and dates), mean and sd
ENVI.Summary <- function(x, Dates) {
  Dates <- as.Date(Dates, format = "%d/%m/%Y")
  min.v <- min(x)
  max.v <- max(x)
  mean.v <- mean(x)
  sd.v <- sd(x)
    min.i <- which.min(x)
    max.i <- which.max(x)
    min.d <- Dates[min.i]
    max.d <- Dates[max.i]
  cat("Min:", min.v,  "on", format(min.d, "%d/%m/%Y"),
      "Max:", max.v, "on", format(max.d, "%d/%m/%Y"),
      "Mean:", mean.v, "SD:", sd.v, "/n")
}

ENVI.Summary(SALINE$Temp, SALINE$Date)
ENVI.Summary(SALINE$DO, SALINE$Date)
ENVI.Summary(SALINE$PH, SALINE$Date)
ENVI.Summary(SALINE$Salinity, SALINE$Date)
ENVI.Summary(SALINE$Turbidity, SALINE$Date)
ENVI.Summary(SALINE$Chlor, SALINE$Date)

################################################################################

# PAIRWISE PERMANOVAS

# PAIRWISE PERMANOVA ON MONTHS
MODEL.M <- pairwise.adonis(T.WATER, SAMPLE$Month)
summary(MODEL.M)
MODEL.M

# PAIRWISE PERMANOVA ON FLOODING PHASES
MODEL.F <- pairwise.adonis(T.WATER, SAMPLE$Flooding_phase)
summary(MODEL.F)
MODEL.F

# PCA ON ENVI
PCA <- prcomp(T.WATER, scale. = TRUE)
# calc percentage
PCA.var <- PCA$sdev^2
PCA.var
PCA.var.per <- round(PCA.var/sum(PCA.var)*100,1)
PCA.var.per

################################################################################

# ANOVAS FOR MONTH AND FLOODING PHASE

# TEMP
Temp.aov.comp <- dredge(aov(Temp ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
Temp.aov <- get.models(Temp.aov.comp, subset = 1)[[1]]
shapiro.test(Temp.aov$residuals)
Temp.aov
summary(Temp.aov)

# DO ~ RESID NOT NORMALLY DISTRIBUTED FOR FIRST MODEL
DO.aov.comp <- dredge(aov(DO ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
DO.aov <- get.models(DO.aov.comp, subset = 1)[[1]]
shapiro.test(DO.aov$residuals)
DO.aov
summary(DO.aov)

# PH
PH.aov.comp <- dredge(aov(PH ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
PH.aov <- get.models(PH.aov.comp, subset = 1)[[1]]
shapiro.test(PH.aov$residuals)
PH.aov
summary(PH.aov)

# SALINITY
Sal.aov.comp <- dredge(aov(Salinity ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
Sal.aov <- get.models(Sal.aov.comp, subset = 1)[[1]]
shapiro.test(Sal.aov$residuals)
Sal.aov
summary(Sal.aov)

# Turbidity
Tur.aov.comp <- dredge(aov(Turbidity ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
Tur.aov <- get.models(Tur.aov.comp, subset = 1)[[1]]
shapiro.test(Tur.aov$residuals)
Tur.aov
summary(Tur.aov)

# Chlorophyll
Chlor.aov.comp <- dredge(aov(Chlor ~ Month + Flooding_phase, data = SALINE, na.action = na.fail))
Chlor.aov <- get.models(Chlor.aov.comp, subset = 1)[[1]]
shapiro.test(Chlor.aov$residuals)
Chlor.aov
summary(Chlor.aov)
# ---


################################################################################

# ENVI ~ MONTH

# PLOT PCA
new.o.1 <- c("February", "March", "April", "May", "July")
Month <- factor(SAMPLE$Month, levels = new.o.1)
# Plot ordination so that points are coloured and shaped according to the groups of interest.
co=c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300", "#CC6633")

par(mar = c(5,5,4,10), xpd = TRUE)
plot(PCA$x, col=co[Month], pch = 16,
     cex=1.7, las = 1, cex.lab = 1.5, xlab = "PC1 (73.9%)", ylab = "PC2 (19.6%)")	
legend('topright', inset = c(-0.2,0), legend = new.o.1, pch = 16, col = co, cex = 1.5, bty = "n")

# Add elipses around groups
ordiellipse(PCA, Month, display = "sites", kind = "sd", col = co, lable = FALSE,
            draw = "polygon", alpha = 30)
ordiellipse(PCA, Month, display = "sites", kind = "sd", col = co, lable = FALSE)

# Arrows
# Create function to find central point of all data points
Centroid <- function(points) {colMeans(points)}

# Apply function to all months
Centroids <- tapply(1:nrow(PCA$x), Month, function(i) Centroid(PCA$x[i, ]))

# Define function for points slightly away from centroid
offset.p <- function(x1, y1, x2, y2, scale = 0.005) {  
  # Calc the direction vector
  dx <- x2 - x1
  dy <- y2 - y1
  # Calc length of direction vector
  len <- sqrt(dx^2 + dy^2)
  # Scale the offset
  if (len == 0) {
    return(c(x1, y1))
  }
  return(c(x1 + scale * dx / len, y1 + scale * dy / len))
}

# Make a gap inbetween arrows (arrows slightly off of centroid point)
# Plot arrows with slight offset
for (i in 1:(length(Centroids) - 1)) {
  start_point <- Centroids[[i]]
  end_point <- Centroids [[i + 1]]
  # DEBUG
  print(paste("start point", start_point))
  print(paste("end point", end_point))
  
  # Calc offset points
  start_offset <- offset.p(start_point[1], start_point[2], end_point[1], end_point[2],scale = 0.03)
  end_offset <- offset.p(end_point[1], end_point[2], start_point[1], start_point[2], scale = 0.03)
  
  arrows(start_offset[1], start_offset[2], end_offset[1], end_offset[2], 
         col = "black", length = 0.1)
}





################################################################################


################################################################################

# BOX PLOTS FOR ENVI ~ MONTH

# DEFINE ORDERS AND COLOURS FOR PLOTS
new.m.o2 <- factor(SALINE$Month, levels = c("February", "March", "April", "May", "July"))
co=c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300")

# PLOT
par(mfrow = c(3,2), oma = c(1,1,1,1), mar = c(4,4,4,4))
boxplot(Temp ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Temperature (C)", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))
boxplot(DO ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Dissolved Oxygen (PPM)", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))
boxplot(PH ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "pH", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))
boxplot(Salinity ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Salinity (PSU)", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))
boxplot(Turbidity ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Turbidity (NTU)", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))
boxplot(Chlor ~ new.m.o2, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Chlorophyll Flouresence (RFU)", col = co, names = c("Feb", "Mar", "Apr", "May", "Jul"))

# BOX PLOTS FOR ENVI ~ FLOODING
order.FL <- factor(SALINE$Flooding_phase, levels = c("BF","DF", "AF", "REC"))
specco <- c("#333333","#FFFFFF","#CCCCCC","#666666")
par(mfrow = c(2,3))
boxplot(Temp ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "Temperature (C)", col = specco)
boxplot(DO ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "Dissolved Oxygen (PPM)", col = specco)
boxplot(PH ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "pH", col = specco)
boxplot(Salinity ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "Salinity (PSU)", col = specco)
boxplot(Turbidity ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "Turbidity (NTU)", col = specco)
boxplot(Chlor ~ order.FL, data=SALINE, las = 1, cex.lab = 1.2, xlab = "Flooding phases", ylab = "Chlorophyll Flouresence (RFU)", col = specco)


################################################################################

# BOX PLOTS FOR ENVI ~ SEASON
order.seas <- factor(SALINE$Season, levels = c("Winter", "Spring", "Summer"))
specco <- c("#006699","#669900","#993300")
par(mfrow = c(2,3))
boxplot(Temp ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "Temperature (C)", col = specco)
boxplot(DO ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "Dissolved Oxygen Concentration", col = specco)
boxplot(PH ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "pH", col = specco)
boxplot(Salinity ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "Salinity (PSU)", col = specco)
boxplot(TDS ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "Total Dissolved Solids", col = specco)
boxplot(Chlor ~ order.seas, data=SALINE, xlab = "Seasons", ylab = "Chlorophyll Concentration", col = specco)


################################################################################

# NOT USED

# ENVI ~ FLOODING PHASES

# PLOT PCA
new.o.2 <- c("BF", "DF", "AF", "REC")
SAMPLE$Flooding_phase <- factor(SAMPLE$Flooding_phase, levels = new.o.2)
# Plot ordination so that points are coloured and shaped according to the groups of interest.
co.2 <- c("#993300", "#006699", "#669900", "#ff9933")

par(mar = c(5,5,4,10), xpd = TRUE)
plot(PCA$x, col=co.2[SAMPLE$Flooding_phase], pch = 16,
     cex=1.7, xlab = "PC1 (61.7%)", ylab = "PC2 (20.8%)")	
legend('topright', inset = c(-0.2,0), legend = new.o.2, pch = 16, col = co.2, cex = 1.5, bty = "n")

# Add elipses around groups
ordiellipse(PCA, SAMPLE$Flooding_phase, display = "sites", kind = "sd", col = co.2, lable = FALSE,
            draw = "polygon", alpha = 30)
ordiellipse(PCA, SAMPLE$Flooding_phase, display = "sites", kind = "sd", col = co.2, lable = FALSE)

# Arrows
# Create function to find central point of all data points
Centroid <- function(points) {colMeans(points)}

# Apply function to all months
Centroids <- tapply(1:nrow(PCA$x), SAMPLE$Flooding_phase, function(i) Centroid(PCA$x[i, ]))

# Define function for points slightly away from centroid
offset.p <- function(x1, y1, x2, y2, scale = 0.005) {  
  # Calc the direction vector
  dx <- x2 - x1
  dy <- y2 - y1
  # Calc length of direction vector
  len <- sqrt(dx^2 + dy^2)
  # Scale the offset
  if (len == 0) {
    return(c(x1, y1))
  }
  return(c(x1 + scale * dx / len, y1 + scale * dy / len))
}

for (i in 1:(length(Centroids) - 1)) {
  start_point <- Centroids[[i]]
  end_point <- Centroids [[i + 1]]
  start_offset <- offset.p(start_point[1], start_point[2], end_point[1], end_point[2],scale = 0.03)
  end_offset <- offset.p(end_point[1], end_point[2], start_point[1], start_point[2], scale = 0.03)
  arrows(start_offset[1], start_offset[2], end_offset[1], end_offset[2], 
         col = "black", length = 0.1)
}


