# Clear R memory
rm(list=ls())
# Close previous graphics windows
graphics.off()
# ---

################################################################################

# LOAD ALL PACKAGES
library(vegan)
library(ggplot2)
library(pairwiseAdonis)
library(ggord)
library(betareg)
library(MuMIn)
library("lme4")

library(lme4)
library(Matrix)
install.packages("lme4", dependencies = TRUE)
library(lme4)
install.packages("Matrix")

install.packages("lmerTest")
library(lmerTest)
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


# PERMANOVA FOR FLOODING PHASES ACROSS ALL SEASONS
PERM.F <- adonis2(T.ABUN ~ Flooding_phase, data=SAMPLE, permutations=999, method="bray")
PERM.F    # Not significant
# ---
PERM.F <- pairwise.adonis(T.ABUN, factors = SAMPLE$Flooding_phase)
?pairwise.adonis

################################################################################

# ANOVAS FOR TOTAL ABUN AND DIV

Total.abun.comp <- dredge(aov(T.TOTAL ~ Month * Flooding_phase, data = SALINE, na.action = na.fail))
Totalabun.mod <- get.models(Total.abun.comp, subset = 1)[[1]]
shapiro.test(Totalabun.mod$residuals)
Totalabun.mod
summary(Totalabun.mod)

TukT <- TukeyHSD(Totalabun.mod)
print(TukT)
summary(TukT)

Total.div.comp <- dredge(aov(TOTAL.shan ~ Month * Flooding_phase, data = SALINE, na.action = na.fail))
Totaldiv.mod <- get.models(Total.div.comp, subset = 2)[[1]]
shapiro.test(Totaldiv.mod$residuals)
Totaldiv.mod
summary(Totaldiv.mod)

################################################################################

# BOX PLOT TOTAL ABUN TOTAL DIV

new.m.o <- c("February", "March", "April", "May", "July")
order.M <- factor(SALINE$Month, levels = c("February", "March", "April", "May", "July"))
Month <- factor(SAMPLE$Month, levels = new.m.o)
co <- c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300")
order.FL <- factor(SALINE$Flooding_phase, levels = c("BF","DF", "AF", "REC"))
specco <- c("#333333","#FFFFFF","#CCCCCC","#666666")

par(mfrow = c(1,2))
boxplot(T.TOTAL ~ order.M, data = SALINE, col = co, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Total Abundance (square root transformed)")
boxplot(TOTAL.shan ~ order.M, data = SALINE, col = co, las = 1, cex.lab = 1.2, xlab = "Month", ylab = "Total Diversity")

par(mfrow= c(1,2))
boxplot(T.TOTAL ~ order.FL, data = SALINE, col = specco, las = 1, cex.lab = 1.2, xlab = "Flooding phase", ylab = "Total abundance (square root transformed)")
boxplot(TOTAL.shan ~ order.FL, data = SALINE, col = specco, las = 1, cex.lab = 1.2, xlab = "Flooding phase", ylab = "Total Diversity")

################################################################################

# BOX PLOTS EACH SPECIES

# SUBSET BY FLOODING PHASE
SALINE$Flooding_phase <- as.factor(SALINE$Flooding_phase)
BF <- SALINE[SALINE$Flooding_phase == "BF",]
DF <- SALINE[SALINE$Flooding_phase == "DF",]
AF <- SALINE[SALINE$Flooding_phase == "AF",]
REC <- SALINE[SALINE$Flooding_phase == "REC",]

BF.ABUN <- BF[,13:29]
BF.M <- colMeans(BF.ABUN)
DF.ABUN <- DF[,13:29]
DF.M <- colMeans(DF.ABUN)
AF.ABUN <- AF[,13:29]
AF.M <- colMeans(AF.ABUN)
REC.ABUN <- REC[,13:29]
REC.M <- colMeans(REC.ABUN)

max.val <- max(T.ABUN)
spec.names <- c("Goby spp.", "Mullet spp.", "Shrimp spp.", "Thin lipped grey mullet", "Amphipod spp.", "Sole", "Golden grey mullet",  "Flatfish spp.", "European flounder", "Cranefly larvae", "Three spined stickleback", "European green crab", "European eel", "European mud scud", "Cockle spp", "blood worm", "polychaete spp.")
specco2 <- c("#660000","#FF0000","#CC3300","#FF9900","#FFFF33","#CCCC00","#669900","#336600","#009966","#33CCCC","#336699","#333399","#330066", "#CC99FF", "#FFCCFF", "#993366", "#990033")
# ---


layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.5))
par(oma = c(5,5,5,5))
par(mar = c(4,4,4,4))
boxplot(BF.ABUN, col = specco2, ylim = c(0, max.val),  las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "Before Flooding")
boxplot(DF.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "During Flooding")
boxplot(AF.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "After Flooding")
boxplot(REC.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "Recovery")
# LEGEND
par(mar = c(0,0,0,0))
plot.new()
legend("center",inset = c(-0.1, 0), legend = spec.names, fill = specco2, title = "Species", cex = 1, ncol = 6)
# ---


# SUBSET BY MONTH
Feb   <- SALINE[SALINE$Month == "February",]
March <- SALINE[SALINE$Month == "March",]
April <- SALINE[SALINE$Month == "April",]
May   <- SALINE[SALINE$Month == "May",]
July  <- SALINE[SALINE$Month == "July",]


Feb.ABUN <- Feb[,13:29]
March.ABUN <- March[,13:29]
April.ABUN <- April[,13:29]
May.ABUN <- May[,13:29]
July.ABUN <- July[,13:29]

par(mfrow = c(3,2))
par(oma = c(1,1,1,1))
par(mar = c(4,4,4,4))
boxplot(Feb.ABUN, col = specco2, ylim = c(0, max.val),  las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "February")
boxplot(March.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "March")
boxplot(April.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "April")
boxplot(May.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "May")
boxplot(July.ABUN, col = specco2, ylim = c(0, max.val), las = 1, cex.lab = 1.2, xlab = "Species", ylab = "Abundance (Helinger transformed)", main = "July")
par(mar = c(0,0,0,0))
plot.new()
legend("center",inset = c(-0.1, 0), legend = spec.names, fill = specco2, title = "Species", cex = 1.4, ncol = 2)

