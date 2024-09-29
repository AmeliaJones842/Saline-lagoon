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

# INVERTEBRATE GENERAL LINEAR MODELS

# COMPARE ALL POSSIBLE MODELS FOR INVERTEBRATE ABUNDANCE
Inv.mod.comp1 <- dredge(glm(T.INV ~ Temp + DO + PH + Salinity + Turbidity + Chlor, data = SALINE, family=gaussian, na.action = na.fail))

Invabun.mod <- get.models(Inv.mod.comp1, subset = 1)[[1]]
shapiro.test(Invabun.mod$residuals)
Invabun.mod
summary(Invabun.mod)

# COMPARE ALL POSSIBLE MODELS FOR INVERTEBRATE DIVERSITY
Inv.mod.comp2 <- dredge(glm(INV.shan ~ Temp + DO + PH + Salinity + Turbidity + Chlor, data = SALINE, family=gaussian, na.action = na.fail))

Invdiv.mod <- get.models(Inv.mod.comp2, subset = 1)[[1]]
shapiro.test(Invdiv.mod$residuals)
Invdiv.mod
summary(Invdiv.mod)

# COMPARE ALL POSSIBLE MODELS FOR SHRIMP ABUNDANCE
Shrimp.mod.comp <- dredge(glm(T.SHRIMP ~ Temp + DO + PH + Salinity + Turbidity + Chlor, data = SALINE, family=gaussian, na.action = na.fail))

Shrimp.mod <- get.models(Shrimp.mod.comp, subset = 1)[[1]]
shapiro.test(Shrimp.mod$residuals)
Shrimp.mod
summary(Shrimp.mod)


################################################################################


# CORRECT SCATTER GRAPHS BASED ON GLM
# INV ABUN/DIV GRAPHS
# GLM STATS TO PRESENT
# RSQUARED FOR ABUNDANCE MODEL
null_model_abun <- glm(T.INV ~ 1, family = gaussian, data = SALINE, na.action = na.fail)
abun_r2 <- 1 - logLik(Invabun.mod) / logLik(null_model_abun)
abun.int <- Inv.mod.comp1[1, "(Intercept)"]
abun.df <- Inv.mod.comp1[1, "df"]

# RSQUARED FOR DIV MODEL
null_model_div <- glm(INV.shan ~ 1, family = gaussian, data = SALINE, na.action = na.fail)
div_r2 <- 1 - logLik(Invdiv.mod) / logLik(null_model_div)
div.int <- Inv.mod.comp2[1, "(Intercept)"]
div.df <- Inv.mod.comp2[1, "df"]

# LEGEND TEXT 1
legend_text1 <- paste0("Inv abun model: glm(INV.abun ~ Temp + DO + PH + Salinity + Turbidity).  R² = ", round(abun_r2, 2), 
                       ",  Intrc ", round(abun.int, 2), 
                       ",  df = ", round(abun.df, 1), "\n")

legend_text2 <- paste0("Inv div model: glm(INV.Shan ~ Chlor).  R² = ", round(div_r2, 2), 
                       ",  Intrc = ", round(div.int, 2), 
                       ",  df = ", round(div.df, 1))

par(mfrow = c(3,2), oma = c(5,0,0,0), mar = c(4,4,2,3), xpd = FALSE)

par(xpd = FALSE)
plot(SALINE$Temp, SALINE$T.INV, pch = 16, las = 1, cex.lab = 1.1, xlab = "Temperature (C)", ylab = "Inv abundance (sq root transformed)", col = "#000000")
Temp.PE <- effect("Temp", Invabun.mod, partial_residuals = TRUE)
lines(Temp.PE$x$Temp, as.vector(Temp.PE$fit), col = "#990000", lwd = 2, lty = 1)

plot(SALINE$DO, SALINE$T.INV, pch = 16, las = 1, cex.lab = 1.1, xlab = "Dissolved Oxygen (PPM) *", ylab = "Inv abundance (sq root transformed)")
DO.PE <- effect("DO", Invabun.mod, partial_residuals = TRUE)
lines(DO.PE$x$DO, as.vector(DO.PE$fit), col = "#669999", lwd = 2, lty = 1)

plot(SALINE$PH, SALINE$T.INV, pch = 16, las = 1, cex.lab = 1.1, xlab = "pH ", ylab = "Inv abundance (sq root transformed)")
PH.PE <- effect("PH", Invabun.mod, partial_residuals = TRUE)
lines(PH.PE$x$PH, as.vector(PH.PE$fit), col = "#CC99FF", lwd = 2, lty = 1)

plot(SALINE$Salinity, SALINE$T.INV, pch = 16, las = 1, cex.lab = 1.1, xlab = "Salinity (PSU) *", ylab = "Inv abundance (sq root transformed)")
Salinity.PE <- effect("Salinity", Invabun.mod, partial_residuals = TRUE)
lines(Salinity.PE$x$Salinity, as.vector(Salinity.PE$fit), col = "#FFCC00", lwd = 2, lty = 1)

plot(SALINE$Turbidity, SALINE$T.INV, pch = 16, las = 1, cex.lab = 1.1, xlab = "Turbidity (NTU)", ylab = "Inv abundance (sq root transformed)")
Turbidity.PE <- effect("Turbidity", Invabun.mod, partial_residuals = TRUE)
lines(Turbidity.PE$x$Turbidity, as.vector(Turbidity.PE$fit), col = "#CC6600", lwd = 2, lty = 1)

plot(SALINE$Chlor, SALINE$INV.shan, pch = 16, las = 1, cex.lab = 1.1, xlab = "Chlorophyll Flouresence (RFU) **", ylab = "Inv Shannon diversity")
Chlor.PE <- effect("Chlor", Invdiv.mod, partial_residuals = TRUE)
lines(Chlor.PE$x$Chlor, as.vector(Chlor.PE$fit), col = "#336600", lwd = 2, lty = 1)

par(xpd = NA) 
legend(x = min(par('usr')[1:2]), y = -0.45, legend = legend_text1, horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.61)
legend(x = min(par('usr')[1:2]), y = -0.6, legend = legend_text2, horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.61)


# SHRIMP ABUN SCATTER 
# LEGEND TEXT 1
null_model_shrimp <- glm(T.SHRIMP ~ 1, family = gaussian, data = SALINE, na.action = na.fail)
shrimp_r2 <- 1 - logLik(Shrimp.mod) / logLik(null_model_shrimp)
shrimp.int <- Shrimp.mod.comp[1, "(Intercept)"]
shrimp.df <- Shrimp.mod.comp[1, "df"]

legend_text1 <- paste0("Shrimp abun model: glm(T.SHRIMP ~ Turbidity + Chlor).  R² = ", round(shrimp_r2, 2), 
                       ",  Intrc ", round(shrimp.int, 2), 
                       ",  df = ", round(shrimp.df, 1), "\n")

par(mfrow = c(1,2), oma = c(5,0,0,0), mar = c(4,4,2,3), xpd = FALSE)
plot(SALINE$Turbidity, SALINE$T.SHRIMP, pch = 16, las = 1, cex.lab = 1.1, xlab = "Turbidity (NTU) *", ylab = "Shrimp Abundance")
Turbidity.PE <- effect("Turbidity", Shrimp.mod, partial_residuals = TRUE)
lines(Turbidity.PE$x$Turbidity, as.vector(Turbidity.PE$fit), col = "#CC6600", lwd = 2, lty = 1)

plot(SALINE$Chlor, SALINE$T.SHRIMP, pch = 16, las = 1, cex.lab = 1.1, xlab = "Chlorophyll Flouresence (RFU) *", ylab = "Shrimp Abundance")
Chlor.PE <- effect("Chlor", Shrimp.mod, partial_residuals = TRUE)
lines(Chlor.PE$x$Chlor, as.vector(Chlor.PE$fit), col = "#336600", lwd = 2, lty = 1)

par(xpd = NA) 
legend(x = min(par('usr')[1:2]), y = -11.9, legend = legend_text1, horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.61)
# ---




# EXTRA STATS


# Correlation of invertebrate and fish
cor.test(SALINE$T.INV, SALINE$T.FISH, method = "kendall")


# ABUN ~ SALINITY
Xvals <- Salinity.PE$x$Salinity
Yvals <-as.vector(Salinity.PE$fit)
COEF <- coef(Invabun.mod)
COEF <- COEF["Salinity"]
intvals <- (Yvals/Xvals) - COEF
Intercept <- mean(intvals)
Lowval <- Intercept + (COEF * 2)
Highval <- Intercept + (COEF * 6)
percchange <- ((Highval - Lowval) / Lowval) * 100
percchange
