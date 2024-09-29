# Clear R memory
rm(list=ls())
# Close previous graphics windows
graphics.off()
# ---

################################################################################

# LOAD ALL PACKAGES
library(vegan)
library(ggplot2)
library(devtools)
library(pairwiseAdonis)
library(ggord)
library(betareg)
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


MFL <- SALINE[,c(2,4)]
WATER2 <- cbind(MFL, WATER.STZ)
WATER.DF <- as.data.frame(WATER2)


# List of independent variables
predictors <- c("Temp", "DO", "PH", "Salinity", "Turbidity", "Chlor")

# Generate all possible combinations of predictors
combinations <- unlist(lapply(1:length(predictors), function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)

# Initialize a dataframe to store model results
model_results <- data.frame(Model = character(), Adjusted_R2 = numeric(), stringsAsFactors = FALSE)

# Function to fit RDA and extract adjusted R2
fit_rda_model <- function(predictor_combination, response_var, data) {
  formula_str <- paste(response_var, "~", paste(predictor_combination, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  model <- rda(formula_obj, data = data)
  adj_r2 <- RsquareAdj(model)$adj.r.squared
  return(list(formula = formula_str, adj_r2 = adj_r2))
}

# Loop through all combinations and fit models
for (i in seq_along(combinations)) {
  result <- fit_rda_model(combinations[[i]], "T.ABUN", WATER.DF)
  model_results[i, ] <- list(Model = result$formula, Adjusted_R2 = result$adj_r2)
}

# View results sorted by Adjusted_R2
model_results <- model_results[order(-model_results$Adjusted_R2), ]

# Print the top model
print(model_results[1, ])

# No perfect Multicolinearity (All VIFs below 10)
RDA <- rda(T.ABUN ~ Temp + DO + PH + Turbidity, data = WATER.DF)
shapiro.test(residuals(RDA))
VIF <- vif.cca(RDA)     

# Compute distance matrix (e.g., Bray-Curtis)
Dist.matrix <- vegdist(T.ABUN, method = "bray")

# Perform db-RDA
db_rda <- capscale(Dist.matrix ~ Temp + DO + PH + Turbidity, data = WATER.DF)

# Summary of the db-RDA
summary(db_rda)

# Plot the db-RDA
plot(db_rda)

# ANOVA ON RDA
# Which RDAs are significant
anova.cca(db_rda, by= "axis")         # Only RDA1 is significant
# Which explanatory variables are statistically significant
anova.cca(db_rda, by= "terms")        # DO, pH and TDS
# ---

new.m.o <- c("February", "March", "April", "May", "July")
Month <- factor(SAMPLE$Month, levels = new.m.o)
CO <- c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300")

ggord(db_rda, grp_in = Month, axes = c("1", "2"), repel = TRUE, 
      arrow = 0.2, cols = c("#0033CC", "#99CCFF", "#336633", "#99CC66", "#993300"), poly = FALSE, 
      polylntyp = "blank")

# NO COLOUR

ggord(db_rda, axes = c("1", "2"), repel = TRUE, 
      arrow = 0.4, cols = ("black"), poly = FALSE, 
      polylntyp = "blank")

################################################################################

# GLMS ON TOTAL ABUNDANCE
# Models
Totalabun.mod1 <- glm(T.TOTAL ~ Temp + DO + PH + Salinity + TDS + Chlor, data = SALINE, family = gaussian())
Totalabun.mod2 <- glm(T.TOTAL ~ Temp + DO  + Salinity + Chlor, data = SALINE, family = gaussian())
Totalabun.mod3 <- glm(T.TOTAL ~ DO + Chlor, data = SALINE, family = gaussian())
Totalabun.mod4 <- glm(T.TOTAL ~ DO * Chlor, data = SALINE, family = gaussian())
Totalabun.mod5 <- glm(T.TOTAL ~ Chlor, data = SALINE, family = gaussian())
# AIC
AIC.1 <- AIC(Totalabun.mod1)
AIC.2 <- AIC(Totalabun.mod2)
AIC.3 <- AIC(Totalabun.mod3)
AIC.4 <- AIC(Totalabun.mod4)
AIC.5 <- AIC(Totalabun.mod5)
# Compare AIC
aic.vals <- c(AIC.1, AIC.2, AIC.3, AIC.4, AIC.5)
print(aic.vals)
# Mod 3 best fit, check residuals and print results
shapiro.test(Totalabun.mod3$residuals)
Totalabun.mod3
summary(Totalabun.mod3)  
# ---

# TOTAL DIVERSITY GLMS * ENVIRONMENTAL VARIABLE
# Models
Totaldiv.mod1 <- glm(TOTAL.shan ~ Temp + DO + PH + Salinity + TDS + Chlor, data = SALINE, family=gaussian())
Totaldiv.mod2 <- glm(TOTAL.shan ~ Temp + DO + PH , data = SALINE, family=gaussian())
Totaldiv.mod3 <- glm(TOTAL.shan ~ Salinity + TDS + Chlor, data = SALINE, family=gaussian())
Totaldiv.mod4 <- glm(TOTAL.shan ~ Temp + DO, data = SALINE, family=gaussian())
Totaldiv.mod5 <- glm(TOTAL.shan ~ Salinity + TDS, data = SALINE, family=gaussian())
Totaldiv.mod6 <- glm(TOTAL.shan ~ Salinity + TDS * Chlor, data = SALINE, family=gaussian())
Totaldiv.mod7 <- glm(TOTAL.shan ~ Chlor, data = SALINE, family=gaussian())
Totaldiv.mod8 <- glm(TOTAL.shan ~ TDS + Chlor, data = SALINE, family=gaussian())
# AIC
AIC.1 <- AIC(Totaldiv.mod1)
AIC.2 <- AIC(Totaldiv.mod2)
AIC.3 <- AIC(Totaldiv.mod3)
AIC.4 <- AIC(Totaldiv.mod4)
AIC.5 <- AIC(Totaldiv.mod5)
AIC.6 <- AIC(Totaldiv.mod6)
AIC.7 <- AIC(Totaldiv.mod7)
AIC.8 <- AIC(Totaldiv.mod8)
# Compare AIC
aic.vals <- c(AIC.1, AIC.2, AIC.3, AIC.4, AIC.5, AIC.6, AIC.7, AIC.8)
print(aic.vals)
# Mod 7 is the best fit, check residuals and print results
shapiro.test(Totaldiv.mod7$residuals)
Totaldiv.mod7
summary(Totaldiv.mod7)  
# ---
