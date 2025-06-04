
############################STATISTICAL ANALYSIS################################################
data <- read.csv("/Users/kelsey/Desktop/EppingAnalysis.csv") #this file contains richness counts, forest cover, distance and NDSI per site
#Determine forest cover scale of effect
# glm: Species richness vs Forest Cover
model_50 <- glm(formula = Richness ~ log(fc50m+1), family = "quasipoisson", data = data)
summary(model_50)
model_100 <- glm(formula = Richness ~ log(fc100m+1), family = "quasipoisson", data = data)
summary(model_100)
model_150 <- glm(formula = Richness ~ log(fc150m+1), family = "quasipoisson", data = data)
summary(model_150)
model_200 <- glm(formula = Richness ~ log(fc200m+1), family = "quasipoisson", data = data)
summary(model_200)
model_250 <- glm(formula = Richness ~ log(fc250m), family = "quasipoisson", data = data)
summary(model_250)
model_300 <- glm(formula = Richness ~ log(fc300m), family = "quasipoisson", data = data)
summary(model_300)
model_350 <- glm(formula = Richness ~ log(fc350m), family = "quasipoisson", data = data)
summary(model_350)
model_400 <- glm(formula = Richness ~ log(fc400m), family = "quasipoisson", data = data)
summary(model_400)
model_450 <- glm(formula = Richness ~ log(fc450m), family = "quasipoisson", data = data)
summary(model_450)
model_500 <- glm(formula = Richness ~ log(fc500m), family = "quasipoisson", data = data)
summary(model_500)
model_550 <- glm(formula = Richness ~ log(fc550m), family = "quasipoisson", data = data)
summary(model_550)
model_600 <- glm(formula = Richness ~ log(fc600m), family = "quasipoisson", data = data)
summary(model_600)
model_650 <- glm(formula = Richness ~ log(fc650m), family = "quasipoisson", data = data)
summary(model_650)
model_700 <- glm(formula = Richness ~ log(fc700m), family = "quasipoisson", data = data)
summary(model_700)
model_750 <- glm(formula = Richness ~ log(fc750m), family = "quasipoisson", data = data)
summary(model_750)
model_800 <- glm(formula = Richness ~ log(fc800m), family = "quasipoisson", data = data)
summary(model_800)

#highest R2 at 600m

ndsi_analysis <- glm(formula = Richness ~ Mean_NDSI, family = "quasipoisson", data = data)
summary(ndsi_analysis)
dist_analysis <- glm(formula = Richness ~ log(Distance_From_M25), family = "quasipoisson", data = data)
summary(dist_analysis)

#Test for correlation
cor.test(log(data$fc600m), data$Mean_NDSI, method = "pearson") 
cor.test(data$Mean_NDSI, log(data$fc600m), method = "pearson") 
cor.test(log(data$fc600m), log(data$Distance_From_M25), method = "pearson") 
cor.test(data$Mean_NDSI, log(data$Distance_From_M25), method = "pearson")

#Plots for report
library(ggplot2)
ggplot(data, aes(x = log(fc600m), y = Richness)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "quasipoisson"),
              se = TRUE, color = "grey40") +
  labs(
    title = "Forest Cover",
    x = "log(Forest Cover 600m) (%)",
    y = "Species Richness"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggplot(data, aes(x = Mean_NDSI, y = Richness)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "quasipoisson"),
              se = TRUE, color = "grey40") +
  labs(
    title = "Mean NDSI",
    x = "Mean NDSI",
    y = "Species Richness"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggplot(data, aes(x = log(Distance_From_M25), y = Richness)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "quasipoisson"),
              se = TRUE, color = "grey40") +
  labs(
    title = "Distance from M25",
    x = "log(Distance from M25) (m)",
    y = "Species Richness"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )


##PLOT SCALE OF EFFECT CURVE##
library(ggplot2)
library(dplyr)

# Load the data
scaleofeffect<- read.csv("new_scale_of_effect.csv") #file of R2 at every forest cover radius
view(scaleofeffect)
dev.off()
# Plot
ggplot(scaleofeffect, aes(x = forest_cover, y = R2_fc)) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "grey40", size = 1.2) +
  geom_vline(xintercept = 600, linetype = "dashed", color = "grey40", size = 1) +  # Add vertical line at 600
  labs(
    x = "Forest Cover Radius (m)",
    y = expression(R^2~"(Species Richness)")
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )


####COMMUNITY COMPOSITION#####
##CC STEP 1 - CREATE PRESENCE/ABSENCE MATRIX
site_species_matrix <- read.csv("/Users/kelsey/Desktop/sitexspeciesmatrix.csv")
sites <- read.csv("EppingAnalysis.csv", header = TRUE)
#Load required libraries
library(dplyr)
library(tidyr)

# Create a clean abundance matrix for community composition
df <- site_species_matrix

# Remove non-species columns
df <- df[, !(names(df) %in% c("X", "Site"))]

# Assign site names as row names
rownames(df) <- site_species_matrix$Site

# Confirm structure (should be all numeric now, no Site column)
str(df)

##CC STEP 2 Calculating dissimilarity indices
#compare 2 sites using total no of species in each site (A+B) and number of species shared (J)
#calculated with vegdist
# Compute Sørensen dissimilarity (Bray with binary=TRUE)
library(vegan)
bray <- vegdist(df, method = "bray", binary = TRUE)
bray_r <- rast(as.matrix(bray))
plot(bray_r, col=hcl.colors(20))

##STEP 3 CC - PERFORM PRINCIPAL COMPONENT ANALYSIS
#takes dissimilarity matrix as a set of independent vectors/ordination axes
pcoa <- cmdscale(bray, k=5, eig=TRUE) 
#graphical CC - 2 dimensions
#Extract the scores
pcoa_axes <- pcoa$points
colnames(pcoa_axes) <- paste0('bray_raw_pcoa_', 1:5)
#check correlations between axes - should not be correlated
zapsmall(cor(pcoa_axes), digits=10) 

#now merge them onto site data to make sure landscape and community metrics match by site
# Convert the pcoa axis values to a data frame and label by site
pcoa_axes_df <- data.frame(pcoa_axes)
pcoa_axes_df$Site <- rownames(pcoa_axes)
# Merge onto the sites data - the PCoA axes get the labels X1 to X5
sites <- merge(sites, pcoa_axes_df, by='Site')
head(sites$Site)
head(pcoa_axes_df$Site)


#STEP 4 - CHECK HOW MANY AXES TO CONSIDER
par(mfrow=c(1,2))
eig <- pcoa$eig[pcoa$eig >0] 
eig
barplot(eig / sum(eig), main='Axis variation')
barplot(cumsum(eig)/ sum(eig), main='Cumulative variation')

# Print the percentage variation of the first 8 
head(sprintf('%0.2f%%', (eig / sum(eig)) * 100), n=8)

##SCALE OF EFFECT CC FOREST COVER
# Model community composition as a function of forest cover
mod_50 <- lm(bray_raw_pcoa_1 ~ log(fc50m+1), data=sites)
summary(mod_50)
# Model community composition as a function of forest cover
mod_100 <- lm(bray_raw_pcoa_1 ~ log(fc100m+1), data=sites)
summary(mod_100)
# Model community composition as a function of forest cover
mod_150 <- lm(bray_raw_pcoa_1 ~ log(fc150m+1), data=sites)
summary(mod_150)
# Model community composition as a function of forest cover
mod_200 <- lm(bray_raw_pcoa_1 ~ log(fc200m+1), data=sites)
summary(mod_200)
# Model community composition as a function of forest cover
mod_250 <- lm(bray_raw_pcoa_1 ~ log(fc250m), data=sites)
summary(mod_250)
# Model community composition as a function of forest cover
mod_300 <- lm(bray_raw_pcoa_1 ~ log(fc300m), data=sites)
summary(mod_300)
# Model community composition as a function of forest cover
mod_350 <- lm(bray_raw_pcoa_1 ~ log(fc350m), data=sites)
summary(mod_350)
# Model community composition as a function of forest cover
mod_400 <- lm(bray_raw_pcoa_1 ~ log(fc400m), data=sites)
summary(mod_400)
# Model community composition as a function of forest cover
mod_450 <- lm(bray_raw_pcoa_1 ~ log(fc450m), data=sites)
summary(mod_450)
# Model community composition as a function of forest cover
mod_500 <- lm(bray_raw_pcoa_1 ~ log(fc500m), data=sites)
summary(mod_500)
# Model community composition as a function of forest cover
mod_550 <- lm(bray_raw_pcoa_1 ~ log(fc550m), data=sites)
summary(mod_550)
# Model community composition as a function of forest cover
mod_600 <- lm(bray_raw_pcoa_1 ~ log(fc600m), data=sites)
summary(mod_600)
# Model community composition as a function of forest cover
mod_650 <- lm(bray_raw_pcoa_1 ~ log(fc650m), data=sites)
summary(mod_650)
# Model community composition as a function of forest cover
mod_700 <- lm(bray_raw_pcoa_1 ~ log(fc700m), data=sites)
summary(mod_700)
# Model community composition as a function of forest cover
mod_750 <- lm(bray_raw_pcoa_1 ~ log(fc750m), data=sites)
summary(mod_750)
# Model community composition as a function of forest cover
mod_800 <- lm(bray_raw_pcoa_1 ~ log(fc800m), data=sites)
summary(mod_800)

#highest R2 at 650m

##PLOT SCALE OF EFFECT CURVE##
library(ggplot2)
library(dplyr)

# Load the data
scaleofeffect<- read.csv("new_scale_of_effect.csv")

ggplot(scaleofeffect, aes(x = forest_cover, y = R2_pcoa)) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "grey40", size = 1.2) +
  geom_vline(xintercept = 650, linetype = "dashed", color = "grey40", size = 1) +  # Add vertical line at 600
  labs(
    x = "Forest Cover Radius (m)",
    y = expression(R^2~"(PCoA 1)")
  ) +
  theme_light() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Plot the model
ggplot(sites, aes(x = log(fc650m), y = bray_raw_pcoa_1)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", se = TRUE) +
  labs(
    title = "Forest Cover",
    x = "log(Forest Cover 650m) (%)",
    y = "PCoA 1 (58.88%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Model community composition as a function of mean NDSI
mod_NDSI <- lm(bray_raw_pcoa_1 ~ Mean_NDSI, data=sites)
summary(mod_NDSI)
mod_m25 <- lm(bray_raw_pcoa_1 ~ log(Distance_From_M25), data=sites)
summary(mod_m25)

#Plot it
ggplot(sites, aes(x = Mean_NDSI, y = bray_raw_pcoa_1)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", se = TRUE) +
  labs(
    title = "Mean NDSI",
    x = "Mean NDSI",
    y = "PCoA 1 (58.88%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(sites, aes(x = log(Distance_From_M25), y = bray_raw_pcoa_1)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", se = TRUE) +
  labs(
    title = "Distance from M25",
    x = "log(Distance from M25) (m)",
    y = "PCoA 1 (58.88%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

########################VARIATION PARTITIONING###############################
#COMMUNITY COMPOSITION PCOA 1
#FIT MODELS & GET adj. R²
# Single sets
mod_1 <- lm(bray_raw_pcoa_1 ~ log(fc650m), data=sites)
summary(mod_1)
#R2 = 0.924
mod_2 <- lm(bray_raw_pcoa_1 ~ Mean_NDSI, data=sites)
summary(mod_2)
#R2 = 0.598
mod_3 <- lm(bray_raw_pcoa_1 ~ log(Distance_From_M25), data=sites)
summary(mod_3)
#R2 = 0.951
# Pairwise combinations
mod_12  <- lm(bray_raw_pcoa_1 ~ log(fc650m) + Mean_NDSI, data=sites)
summary(mod_12)
#R2 = 0.923
mod_13  <- lm(bray_raw_pcoa_1 ~ log(fc650m) + log(Distance_From_M25), data=sites)
summary(mod_13)
#R2 = 0.978
mod_23  <- lm(bray_raw_pcoa_1 ~ Mean_NDSI   + log(Distance_From_M25), data=sites)
summary(mod_23)
#R2 = 0.945
# Full model (all three)
mod_123 <- lm(bray_raw_pcoa_1 ~ log(fc650m) + Mean_NDSI + log(Distance_From_M25), data = sites)
summary(mod_123)
#R2 = 0.975

R2_1   <- 0.924
R2_2   <- 0.598
R2_3   <- 0.951
R2_12  <- 0.923
R2_13  <- 0.978
R2_23  <- 0.945
R2_123 <- 0.975

#Calculate partitions
# Shared three-way overlap
g <- R2_1 + R2_2 + R2_3 - R2_12 - R2_13 - R2_23 + R2_123

# Shared two-way overlaps
d <- R2_12 - R2_1 - R2_2 + g        # X1 ∩ X2 only
e <- R2_13 - R2_1 - R2_3 + g          # X1 ∩ X3 only
f <- R2_23 - R2_2 - R2_3 + g          # X2 ∩ X3 only

# Unique fractions
a <- R2_123 - R2_23  # unique to X1 (fc650m)
b <- R2_123 - R2_13  # unique to X2 (NDSI)
c <- R2_123 - R2_12  # unique to X3 (Dist from M25)
# Residual (unexplained)
h <- 1 - R2_123

# ---- 3. DISPLAY RESULT ----
parts <- round(c(
  unique_forest_cover = a,
  unique_mean_NDSI    = b,
  unique_dist_M25     = c,
  shared_forest_NDSI  = -d,
  shared_forest_dist  = -e,
  shared_NDSI_dist    = -f,
  shared_all_three    = g,
  unexplained         = h
), 4)

print(parts)


#RICHNESS
data <- read.csv("/Users/kelsey/Desktop/EppingAnalysis.csv")

#Pseudo-R2 from GLMs
# Single sets
mod_1 <- glm(formula = Richness ~ log(fc600m), family = "quasipoisson", data = data)
summary(mod_1)
#R2 = 0.892
mod_2 <- glm(formula = Richness ~ Mean_NDSI, family = "quasipoisson", data = data)
summary(mod_2)
#R2 = 0.402
mod_3 <- glm(formula = Richness ~ log(Distance_From_M25), family = "quasipoisson", data = data)
summary(mod_3)
#R2 = 0.864
# Pairwise combinations
mod_12  <- glm(formula = Richness ~ log(fc600m) + Mean_NDSI, family = "quasipoisson", data = data)
summary(mod_12)
#R2 = 0.895
mod_13  <-  glm(formula = Richness ~ log(fc600m) + log(Distance_From_M25), family = "quasipoisson", data = data)
summary(mod_13)
#R2 = 0.927
mod_23  <-  glm(formula = Richness ~ Mean_NDSI + log(Distance_From_M25), family = "quasipoisson", data = data)
summary(mod_23)
#R2 = 0.893

# Full model (all three)
mod_123 <- glm(formula = Richness ~ log(fc600m) + Mean_NDSI + log(Distance_From_M25), family = "quasipoisson", data = data)
summary(mod_123)
#R2 = 0.952
# Extract adjusted R²
getR2 <- function(m) summary(m)$adj.r.squared
R2_1   <- 0.892
R2_2   <- 0.402
R2_3   <- 0.864
R2_12  <- 0.895
R2_13  <- 0.927
R2_23  <- 0.893
R2_123 <- 0.952


# Shared three-way overlap
g <- R2_1 + R2_2 + R2_3 - R2_12 - R2_13 - R2_23 + R2_123

# Shared two-way overlaps
d <- R2_12 - R2_1 - R2_2 + g        # X1 ∩ X2 only
e <- R2_13 - R2_1 - R2_3 + g          # X1 ∩ X3 only
f <- R2_23 - R2_2 - R2_3 + g          # X2 ∩ X3 only

# Unique fractions
a <- R2_123 - R2_23  # unique to X1 (fc650m)
b <- R2_123 - R2_13  # unique to X2 (NDSI)
c <- R2_123 - R2_12  # unique to X3 (Dist from M25)
# Residual (unexplained)
h <- 1 - R2_123

# ---- 3. DISPLAY RESULT ----
parts <- round(c(
  unique_forest_cover = a,
  unique_mean_NDSI    = b,
  unique_dist_M25     = c,
  shared_forest_NDSI  = -d,
  shared_forest_dist  = -e,
  shared_NDSI_dist    = -f,
  shared_all_three    = g,
  unexplained         = h
), 4)

print(parts)

