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
library(effects)
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

# FISH GENERAL LINEAR MODELS

# COMPARE ALL POSSIBLE MODELS FOR FISH ABUNDANCE
Fish.mod.comp1 <- dredge(glm(T.FISH ~ Temp + DO + PH + Salinity + Turbidity + Chlor, data = SALINE, family=gaussian, na.action = na.fail))

Fishabun.mod <- get.models(Fish.mod.comp1, subset = 1)[[1]]
shapiro.test(Fishabun.mod$residuals)
Fishabun.mod
summary(Fishabun.mod)

# COMPARE ALL POSSIBLE MODELS FOR FISH DIVERSITY
Fish.mod.comp2 <- dredge(glm(FISH.shan ~ Temp + DO + PH + Salinity + Turbidity + Chlor, data = SALINE, family=gaussian, na.action = na.fail))

Fishdiv.mod <- get.models(Fish.mod.comp2, subset = 1)[[1]]
shapiro.test(Fishdiv.mod$residuals)
Fishdiv.mod
summary(Fishdiv.mod)


# FISH ABUNDANCE AND DIVERSITY RESULTS
Fishabun.mod
summary(Fishabun.mod)
Fishdiv.mod
summary(Fishdiv.mod)
# ---


################################################################################

# Fish abundance and diversity SIGNIF scatter plots combined

# GLM STATS TO PRESENT
# RSQUARED FOR ABUNDANCE MODEL
null_model_abun <- glm(T.FISH ~ 1, family = gaussian, data = SALINE, na.action = na.fail)
abun_r2 <- 1 - logLik(Fishabun.mod) / logLik(null_model_abun)
abun.int <- Fish.mod.comp1[1, "(Intercept)"]
abun.df <- Fish.mod.comp1[1, "df"]

# RSQUARED FOR DIVERSITY MODEL
null_model_div <- glm(FISH.shan ~ 1, family = gaussian, data = SALINE, na.action = na.fail)
div_r2 <- 1 - logLik(Fishdiv.mod) / logLik(null_model_div)
div.int <- Fish.mod.comp2[1, "(Intercept)"]
div.df <- Fish.mod.comp2[1, "df"]

# LEGEND TEXT 1
legend_text1 <- paste0("Fish abun model: R² = ", round(abun_r2, 2), 
                       ",  Intrc ", round(abun.int, 2), 
                       ",  df = ", round(abun.df, 1), "\n")

legend_text2 <- paste0("Fish div model: R² = ", round(div_r2, 2), 
                       ",  Intrc = ", round(div.int, 2), 
                       ",  df = ", round(div.df, 1))


par(mfrow = c(2,2), oma = c(5,0,0,0), mar = c(4,4,2,5), xpd = FALSE)

par(xpd = FALSE)
plot(SALINE$Temp, SALINE$T.FISH, pch = 16, las = 1, cex.lab = 1.1, xlab = "Temperature (C)", ylab = "Fish abundance (sq root transformed) **", col = "#000000")
Temp.PE1 <- effect("Temp", Fishabun.mod, partial_residuals = TRUE)
lines(Temp.PE1$x$Temp, as.vector(Temp.PE1$fit), col = "#990000", lwd = 2, lty = 1)
par(new = TRUE)
plot(SALINE$Temp, SALINE$FISH.shan, pch = 1, las = 1, xlab = "", ylab = "", col = "#000000", axes = FALSE)
axis(side = 4, las = 1)
mtext("Fish diversity (Shannons index) ***", side = 4, line = 3, cex.lab = 1.1)
Temp.PE2 <- effect("Temp", Fishdiv.mod, partial_residuals = TRUE)
lines(Temp.PE2$x$Temp, as.vector(Temp.PE2$fit), col = "#990000", lwd = 2, lty = 2)


plot(SALINE$DO, SALINE$T.FISH, pch = 16, las = 1, cex.lab = 1.1, xlab = "Dissolved Oxygen (PPM)", ylab = "Fish abundance (sq root transformed) ***")
DO.PE1 <- effect("DO", Fishabun.mod, partial_residuals = TRUE)
lines(DO.PE1$x$DO, as.vector(DO.PE1$fit), col = "#669999", lwd = 2, lty = 1)
par(new = TRUE)
plot(SALINE$DO, SALINE$FISH.shan, pch = 1, las = 1, xlab = "", ylab = "", col = "#000000", axes = FALSE)
axis(side = 4, las = 1)
mtext("Fish diversity (Shannons index) **", side = 4, line = 3, cex.lab = 1.1)
DO.PE2 <- effect("DO", Fishdiv.mod, partial_residuals = TRUE)
lines(DO.PE2$x$DO, as.vector(DO.PE2$fit), col = "#669999", lwd = 2, lty = 2)

plot(SALINE$Turbidity, SALINE$T.FISH, pch = 16, las = 1, cex.lab = 1.1, xlab = "Turbidity (NTU)", ylab = "Fish abundance (sq root transformed) **")
Turbidity.PE1 <- effect("Turbidity", Fishabun.mod, partial_residuals = TRUE)
lines(Turbidity.PE1$x$Turbidity, as.vector(Turbidity.PE1$fit), col = "#CC6600", lwd = 2, lty = 1)
par(new = TRUE)
plot(SALINE$Turbidity, SALINE$FISH.shan, pch = 1, las = 1, xlab = "", ylab = "", col = "#000000", axes = FALSE)
axis(side = 4, las = 1)
mtext("Fish diversity (Shannons index) **", side = 4, line = 3, cex.lab = 1.1)
Turbidity.PE2 <- effect("Turbidity", Fishdiv.mod, partial_residuals = TRUE)
lines(Turbidity.PE2$x$Turbidity, as.vector(Turbidity.PE2$fit), col = "#CC6600", lwd = 2, lty = 2)

par(xpd = NA)  # Allow legend to be drawn outside plot region
legend(x = min(par('usr')[1:2]), y = -0.4, legend = c("Fish abun partial effect", "Fish abun samples"),
       lty = c(1, 0), pch = c(NA,16), col = c("#000000", "#000000"), horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.1)
legend(x = min(par('usr')[1:2]), y = -0.5, legend = legend_text1, horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.065)

par(xpd = FALSE)
plot(SALINE$Chlor, SALINE$T.FISH, pch = 16, las = 1, cex.lab = 1.1, xlab = "Chlorophyll Flourescence (RFU)", ylab = "Fish abundance (sq root transformed) **")
Chlor.PE1 <- effect("Chlor", Fishabun.mod, partial.residuals = TRUE)
lines(Chlor.PE1$x$Chlor, as.vector(Chlor.PE1$fit), col = "#336600", lwd = 2, lty = 1)
par(new = TRUE)
plot(SALINE$Chlor, SALINE$FISH.shan, pch = 1, las = 1, xlab = "", ylab = "", axes = FALSE)
axis(side = 4, las = 1)
mtext("Fish diversity (Shannons index) *", side = 4, line = 3, cex.lab = 1.1)
Chlor.PE2 <- effect("Chlor", Fishdiv.mod, partial.residuals = TRUE)
lines(Chlor.PE2$x$Chlor, as.vector(Chlor.PE2$fit), col = "#336600", lwd = 2, lty = 2)

par(xpd = NA)  # Allow legend to be drawn outside plot region
legend(x = mean(par('usr')[1:2]), y = -0.4, legend = c("Fish div partial effect", "Fish div samples"),
       lty = c(2, 0), pch = c(NA,1), col = c("#000000", "#000000"), horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.44)
legend(x = min(par('usr')[1:2]), y = -0.5, legend = legend_text2, horiz = TRUE, bty = "n", cex = 1.2, xjust = 0.01)


################################################################################


# EXTRA STATS

# ABUN ~ DO
Xvals <- DO.PE1$x$DO
Yvals <-as.vector(DO.PE1$fit)
COEF <- coef(Fishabun.mod)
COEF <- COEF["DO"]
intvals <- (Yvals/Xvals) - COEF
Intercept <- mean(intvals)
Lowval <- Intercept + (COEF * 3)
Highval <- Intercept + (COEF * 8)
percchange <- ((Highval - Lowval) / Lowval) * 100
percchange

# DIV ~ TEMP
Xvals <- Temp.PE2$x$Temp
Yvals <-as.vector(Temp.PE2$fit)
COEF <- coef(Fishdiv.mod)
COEF <- COEF["Temp"]
intvals <- (Yvals/Xvals) - COEF
Intercept <- mean(intvals)
Lowval <- Intercept + (COEF * 10)
Highval <- Intercept + (COEF * 15)
percchange <- ((Highval - Lowval) / Lowval) * 100
percchange

# DIV ~ Chlor
Xvals <- Chlor.PE2$x$Chlor
Yvals <-as.vector(Chlor.PE2$fit)
COEF <- coef(Fishdiv.mod)
COEF <- COEF["Chlor"]
intvals <- (Yvals/Xvals) - COEF
Intercept <- mean(intvals)
Lowval <- Intercept + (COEF * 0.5)
Highval <- Intercept + (COEF * 2.5)
percchange <- ((Highval - Lowval) / Lowval) * 100
percchange



cor.test(SALINE$FISH, SALINE$FISH.shan)
plot(SALINE$FISH, SALINE$FISH.shan, pch = 16)
fit <- lm(SALINE$FISH.shan ~SALINE$FISH)
abline(fit, col = "red")

