# Clear R memory
rm(list=ls())
# Close previous graphics windows
graphics.off()
# ---

################################################################################

# DATA LOADING

# LOAD ALL PACKAGES
library("vegan")
library("ggplot2")
library("pairwiseAdonis")
library("ggord")
library("betareg")
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

# SEASONALITY

# NON METRIC MULTIDIMENSIONAL SCALING (NMDS)
# Run NMDS to visualise the community structure of each data point (3 seines) uing Bray-Curtis distance as data is count data. 
NMDS <- metaMDS(T.ABUN, distance = "bray", k = 2)
plot(NMDS)
# STRESS PLOT
# Shows us how well the NMDS model fits the data.
SPNMDS <- stressplot(NMDS)
NMDS$stress    
# ---

fit <- envfit(NMDS, T.ABUN, perm = 999)
plot(NMDS)
plot(fit, p.max = 0.05, col = "black")


# PERMANOVA FOR MONTHS
PERM.M1 <- adonis2(T.ABUN ~ Month, data=SAMPLE, permutations=999, method="bray")
PERM.M1    # Significant
# Pairwise to test which groups are significant.
MODEL.M2 <- pairwise.adonis(T.ABUN, SAMPLE$Month)
summary(MODEL.M2)
MODEL.M2
# ---  

new.m.o <- c("February", "March", "April", "May", "July", "August")
Month <- factor(SAMPLE$Month, levels = new.m.o)


CO <- c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300", "#CC6633")

par(mar = c(5,5,4,10), xpd = TRUE)
plot(NMDS$points, col=CO[Month], pch = 16, cex = 1.5, las = 1, cex.lab = 1.5, xlab = "NMDS1", ylab = "NMDS2")
text(x = max(NMDS$points[,1]), y = min(NMDS$points[,2]), 
    labels = paste("Stress =", round(NMDS$stress, 3)), pos = 2, cex = 1.2, col = "black")
    legend("bottomright", inset = c(-0.24,0.55),
           legend = c("S.spp: Shrimp species",
                      "M.spp: Mullet species",
                      "G.spp: Goby species",
                      "EMS: European mud scud",
                      "Am: Amphipod",
                      "EF: European flounder",
                      "EGC: European green crab"), 
                      cex = 0.9, col = "black", bty = "n")
legend('topright', inset = c(-0.2,0), legend = new.m.o, pch = 16, col = CO, cex = 1.5, bty = "n")

# Add circles around each month
ordiellipse(NMDS, Month, display = "sites", kind = "sd", col = CO, lable = FALSE,
            draw = "polygon", alpha = 30)
ordiellipse(NMDS, Month, display = "sites", kind = "sd", col = CO, lable = FALSE)

# Create function to find central point of all data points
Centroid <- function(points) {colMeans(points)}

# Apply function to all months
Centroids <- tapply(1:nrow(NMDS$points), Month, function(i) Centroid(NMDS$points[i, ]))

# Mark centroid of each month with an x
#for (i in 1:length(Centroids)) {
#  points(Centroids[[i]][1], Centroids[[i]][2], pch = 4, col = CO[i], cex = 1)
#}

# Add arrows between each centroid from one month to the next
#for (i in 1:(length(Centroids) - 1)) {
#  arrows(Centroids[[i]][1], Centroids[[i]][2], Centroids[[i + 1]][1], Centroids[[i + 1]][2],
#         col = "black", length = 0.1)
#}

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
  start_offset <- offset.p(start_point[1], start_point[2], end_point[1], end_point[2],scale = 0.005)
  end_offset <- offset.p(end_point[1], end_point[2], start_point[1], start_point[2], scale = 0.005)
  
  arrows(start_offset[1], start_offset[2], end_offset[1], end_offset[2], 
         col = c("black"), length = 0.08)
}

par(mfrow= c(1,1))
plot(fit, p.max = 0.05, col = "#999999", cex = 1.2, lwd = 2)

# ---

