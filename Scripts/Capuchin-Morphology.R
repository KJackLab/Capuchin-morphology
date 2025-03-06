# -- Dominance rank, facial morphology, and testes size in male white-faced capuchins: evidence for pre- and post-mating competition -- #
# Load Required Packages ####
library(lmerTest)
library(lme4)
library(ggplot2)
library(dplyr)
library(DHARMa)
library(tidyverse)
library(car)
library(psych)
library(ggfortify)  # For PCA plotting
library(effects)
library(MuMIn)
library(stats)
library(graphics)
library(zoo)
library(tidyr)

#Load and prepare data for scrotum analysis####
df <- photogrammetry_dataset_stable

# Remove rows with all NA values
df <- df %>% filter(rowSums(is.na(.)) != ncol(.))

# Convert relevant columns to factors and numeric types
df <- df %>%
  mutate(
    ID = as.factor(individual),
    group = as.factor(group),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank)) %>%
  filter(dominance_rank %in% c("Alpha", "Subordinate")) %>%
  mutate(across(c(brow_width, muzzle_width, facial_width, facial_length, scrotum_width, left_test_length, right_test_length, body_length), as.numeric))

# Check  unique levels of  dominance_rank
print(unique(df$dominance_rank))

df$date <- as.Date(df$date, format = "%Y.%m.%d")

df$laser_type <- ifelse(df$date < as.Date("2023-07-01"), "1", "2")


# Define the function to impute missing values
impute_closest <- function(x) {
  non_na_indices <- which(!is.na(x)) # Get indices of non-NA values
  if (length(non_na_indices) == 0) return(x) # Return unchanged if all values are NA
  
  imputed <- sapply(seq_along(x), function(i) {
    if (!is.na(x[i])) {
      return(x[i]) # Keep non-missing values
    } else {
      # Find the closest non-missing value
      distances <- abs(non_na_indices - i)
      closest_index <- non_na_indices[which.min(distances)]
      return(x[closest_index])
    }
  })
  
  return(imputed)
}

# Impute missing body lengths
df <- df %>%
  group_by(ID) %>%
  arrange(date) %>% # Ensure data is sorted by ID and date
  mutate(body_length = impute_closest(body_length)) %>%
  ungroup()

# Check the data structure 
str(df)
head(df)

# Summary statistics 
summary(df)

# Make sure each individual is either alpha or subordinate. NOT both
table(df$individual, df$dominance_rank)

#Descriptive statistics####
# Summary statistics for measurements by dominance rank
df_long <- df %>%
  pivot_longer(cols = c(scrotum_width, facial_width, muzzle_width, brow_width, facial_length),
               names_to = "body_part",
               values_to = "measurement")

summary_stats <- df_long %>%
  group_by(body_part, dominance_rank) %>%
  summarise(
    mean_value = mean(measurement, na.rm = TRUE),
    sd_value = sd(measurement, na.rm = TRUE),
    min_value = min(measurement, na.rm = TRUE),
    max_value = max(measurement, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

print(summary_stats)
 
#Scrotum Width analysis####
# Data cleaning: Remove rows with any NAs in the relevant columns
df_scrotum_width <- df %>%
  filter(!is.na(scrotum_width) & !is.na(dominance_rank) & !is.na(age)& !is.na(Temp) & !is.na(individual))

# count sample size for this analysis
table(df_scrotum_width$individual, df_scrotum_width$dominance_rank)

# Fit mixed-effects model for scrotum width
scrotum_model <- lme4::lmer(scrotum_width ~ dominance_rank + body_length + age + (1 | individual), data = df_scrotum_width)

# Check structure and distribution
hist(df_scrotum_width$scrotum_width)
shapiro.test(df_scrotum_width$scrotum_width)  

# Diagnostics and model output
summary(scrotum_model)
drop1(scrotum_model, test ="Chisq")
r.squaredGLMM(scrotum_model)

# Simulate and plot residuals
res_scrotum_model1 <- simulateResiduals(scrotum_model1)
plot(res_scrotum_model1)  # Residual diagnostics #LOOKS GOOD
plot(allEffects(scrotum_model1))

#test outliers
testOutliers(res_scrotum_model)

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_scrotum <- testUniformity(res_scrotum_model)
# Print KS test result
print(ks_test_result_scrotum)

# Perform dispersion test
dispersion_test_result_scrotum <- testDispersion(res_scrotum_model)
# Print dispersion test result
print(dispersion_test_result_scrotum) #p=0.88

# Plot scrotum width by dominance rank
p1<-ggplot(df_clean_scrotum_width, aes(x = dominance_rank, y = scrotum_width, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Scrotum Width (mm)") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend
p1

# Scatter plot showing scrotum width vs age by dominance rank
a1 <- ggplot(df, aes(x = age, y = scrotum_width, color = dominance_rank, shape = dominance_rank)) +  # Correctly close the aes function
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.7) +  # Adds a single regression line for all data
  labs(
    title = "Scrotum Width (mm)",  # Corrected the missing quotation mark
    x = "Age (years)",
    y = "Scrotum Width (mm)",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensures the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  # Adjusts the legend position to the bottom
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and pale green for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a1

#Load/organize/impute data for Facial analyses####
dff <- photogrammetry_dataset_stable #dataset WITHOUT outlier (LE)

# Remove rows with all NA values
dff <- dff %>% filter(rowSums(is.na(.)) != ncol(.))

# Remove individual 'LE'
df_filtered <- dff %>%
  filter(individual != "LE")

# Convert relevant columns to factors and numeric, and transform data
# Only stable periods.
df_filtered <- df_filtered %>%
  mutate(
    ID = as.factor(individual),
    group = as.factor(group),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank)) %>%
  filter(dominance_rank %in% c("Alpha", "Subordinate")) %>%
  mutate(across(c(brow_width, muzzle_width, facial_width, facial_length, scrotum_width, left_test_length, right_test_length, body_length), as.numeric))

#imputate body length as Last Observation Carried Forward (LOCF)
df_filtered <- df_filtered  %>%
  group_by(ID) %>%
  arrange(date) %>% # Ensure the data is ordered by date
  mutate(
    body_length = if (all(is.na(body_length))) NA else na.locf(body_length, na.rm = FALSE), # Forward fill
    body_length = if (all(is.na(body_length))) NA else na.locf(body_length, fromLast = TRUE) # Backward fill
  ) %>%
  ungroup()

#Check the unique levels of the dominance_rank column
print(unique(df_filtered$dominance_rank))



# Check the structure of the dataset
str(df_filtered)
head(df_filtered)

# Summary statistics of the dataset
summary(df_filtered)

df_filtered$date <- as.Date(df_filtered$date, format = "%Y.%m.%d")

df_filtered$laser_type <- ifelse(df_filtered$date < as.Date("2023-07-01"), "1", "2")

# Make sure each ind. is either alpha or subordinate. NOT both
table(df_filtered$individual, df_filtered$dominance_rank)

# Convert to long format for easier summarization
df_long <- df_filtered %>%
  pivot_longer(cols = c(scrotum_width, facial_width, muzzle_width, brow_width, facial_length),
               names_to = "body_part",
               values_to = "measurement")

# Calculate summary statistics for each body part and dominance rank
summary_stats <- df_long %>%
  group_by(body_part, dominance_rank) %>%
  summarise(
    mean_value = mean(measurement, na.rm = TRUE),
    sd_value = sd(measurement, na.rm = TRUE),
    min_value = min(measurement, na.rm = TRUE),
    max_value = max(measurement, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()
print(summary_stats)

#PCA test for FACIAL models####
# Standardize the data
facedf <- df_filtered[,c("individual","group", "dominance_rank", "date", "brow_width", "muzzle_width", "facial_width", "facial_length")]

# Remove rows with NAs in the facial measurements
facedf_clean <- facedf %>%
  filter(!is.na(brow_width) & !is.na(muzzle_width) & !is.na(facial_width) & !is.na(facial_length))

# Standardize the data
columns_to_standardize <- facedf_clean[, c("brow_width", "muzzle_width", "facial_width", "facial_length")]
standardized_data <- scale(columns_to_standardize)

#print pca results
pca_result <- principal(as.matrix(standardized_data), nfactors=4, rotate="varimax")
print(pca_result)

KMO(standardized_data)
cortest.bartlett(standardized_data)

#Facial width analysis####
# Data Cleaning
df_facial_width <- df_filtered %>%
  filter(!is.na(facial_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_facial_width <- df_facial_width %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_width))

#Check Assumptions with Diagnostics
hist(df_facial_width$facial_width) #looks v normal
shapiro.test(df_facial_width$facial_width) # normal p=0.80

# count sample size for this analysis
table(df_facial_width$individual, df_facial_width$dominance_rank)

# Fit mixed-effects model for facial length
facial_width_model <- lme4::lmer(facial_width ~ dominance_rank + body_length + age + (1 | individual), data = df_facial_width)

#Diagnostics and model output
summary(facial_width_model)
drop1(facial_width_model, test ="Chisq")
r.squaredGLMM(facial_width_model)

# Simulate residuals and plot diagnostics
res_facial_width_model <- simulateResiduals(facial_width_model)
plot(res_facial_width_model)  # looks fine
plot(allEffects(facial_width_model)) #looks fine

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_fw <- testUniformity(res_facial_width_model)
# Print KS test result
print(ks_test_result_fw) 

# Perform dispersion test
dispersion_test_result_fw <- testDispersion(res_facial_width_model)
# Print dispersion test result
print(dispersion_test_result_fw) #fine p=0.63

#test outliers
testOutliers(res_facial_width_model)

# Plot for Facial Width
p2 <-ggplot(df_facial_width, aes(x = dominance_rank, y = facial_width, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Facial Width (mm)") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend

p2

# Scatter plot showing scrotum width vs age by dominance rank
a2<-ggplot(df_facial_width, aes(x = age, y = facial_width, color = dominance_rank, shape = dominance_rank)) +
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.7) +  # Adds a single regression line for all data
  labs(
    title = "Facial Width (mm)",  # Corrected the missing quotation mark
    x = "Age (years)",
    y = "Facial Width (mm)",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensures the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  # Adjusts the legend position to the bottom
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and pale green for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a2

#Muzzle width analysis####
#Data Cleaning
df_muzzle <- df_filtered %>%
  filter(!is.na(muzzle_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_muzzle<- df_muzzle %>%
  filter(complete.cases(dominance_rank, body_length, age, muzzle_width))

#Check Assumptions with Diagnostics
hist(df_muzzle$muzzle_width)
shapiro.test(df_muzzle$muzzle_width) 

# count sample size for this analysis
table(df_muzzle$individual, df_muzzle$dominance_rank)

#Fix mixed-effects model
muzzle_width_model <- lme4::lmer(muzzle_width ~ dominance_rank + body_length + age + (1 | individual), data = df_muzzle)

#Diagnostics and model output
summary(muzzle_width_model)
drop1(muzzle_width_model, test ="Chisq")
r.squaredGLMM(muzzle_width_model)

#Model Diagnostics
res_muzzle_width_model <- simulateResiduals(muzzle_width_model)
plot(res_muzzle_width_model)  
plot(allEffects(muzzle_width_model))

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_muzzle <- testUniformity(res_muzzle_width_model)
# Print KS test result
print(ks_test_result_muzzle)

# Perform dispersion test
dispersion_test_result_muzzle <- testDispersion(res_muzzle_width_model)
# Print dispersion test result
print(dispersion_test_result_muzzle) 

#test outliers
testOutliers(res_muzzle_width_model)

#Plot for Muzzle Width
p3<-ggplot(df_muzzle, aes(x = dominance_rank, y = muzzle_width, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Muzzle Width (mm)") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend
p3

# Scatter plot showing scrotum width vs age by dominance rank
a3<-ggplot(df_muzzle, aes(x = age, y = muzzle_width, color = dominance_rank, shape = dominance_rank)) +
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.7) +  # Adds a single regression line for all data
  labs(
    title = "Muzzle Width (mm)",  # Corrected the missing quotation mark
    x = "Age (years)",
    y = "Muzzle Width (mm)",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensures the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  # Adjusts the legend position to the bottom
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and pale green for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a3

#Brow width analysis####
# Filter out missing values for brow_width
df_brow <- df_filtered %>%
  filter(!is.na(brow_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_brow <- df_brow %>%
  filter(complete.cases(dominance_rank, body_length, age, brow_width))

#Check Assumptions 
hist(df_brow$brow_width)
shapiro.test(df_brow$brow_width) 

# count sample size for this analysis
table(df_brow$individual, df_brow$dominance_rank)

#Fit mixed-effects model
brow_model <- lme4::lmer(brow_width ~ dominance_rank + body_length + age + (1 | individual), data = df_brow)

#Diagnostics and model output
summary(brow_model)
drop1(brow_model, test="Chisq")
r.squaredGLMM(brow_model)

#Residuals diagnostics
residuals_sim <- simulateResiduals(fittedModel = brow_model)
plot(residuals_sim)

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_brow <- testUniformity(res_brow_width_model)
# Print KS test result
print(ks_test_result_brow) #p=0.53

# Perform dispersion test
dispersion_test_result_brow <- testDispersion(res_brow_width_model)
# Print dispersion test result
print(dispersion_test_result_brow) #p=0.632

#test outliers
testOutliers(res_brow_width_model)

# Plot for Brow Width by dominance rank
p4<-ggplot(df_brow, aes(x = dominance_rank, y = brow_width, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Brow Width (mm)") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend
p4

# Scatter plot showing Brow Width vs age by dominance rank
a4 <- ggplot(df_brow, aes(x = age, y = brow_width, color = dominance_rank, shape = dominance_rank)) +
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.5) +  # Adds a single regression line for all data
  labs(
    title = "Brow Width (mm)",
    x = "Age (years)",
    y = "Brow Width (mm)",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensure the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and gold for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a4

#Facial Height analysis#### 
# Remove rows with missing values for the columns used in the model
df_length <- df_filtered %>%
  filter(!is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_length <- df_length %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_length))

# Check assumptions
hist(df_length$facial_length)
shapiro.test(df_length$facial_length)  # normal (p=.71)

# count sample size for this analysis
table(df_length$individual, df_length$dominance_rank)

# Fit mixed-effects model for facial length
facial_length_model <- lme4::lmer(facial_length ~ dominance_rank + body_length + age + (1 | individual), data = df_length)

#Diagnostics and model output
summary(facial_length_model)
drop1(facial_length_model, test= "Chisq")
r.squaredGLMM(facial_length_model)

# Simulate residuals and plot diagnostics
res_facial_length_model <- simulateResiduals(facial_length_model)
plot(res_facial_length_model)  # Residual diagnostics
plot(allEffects(facial_length_model))

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_length <- testUniformity(res_facial_length_model)
print(ks_test_result_length) 

# Perform dispersion test
dispersion_test_result_length <- testDispersion(res_facial_length_model)
# Print dispersion test result
print(dispersion_test_result_length) 

#test outliers
testOutliers(res_facial_length_model)

# Plot for facial length by dominance rank
p5<-ggplot(df_length, aes(x = dominance_rank, y = facial_length, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Facial Height (mm)") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend
p5

# Scatter plot showing facial length vs age by dominance rank
a5 <- ggplot(df_length, aes(x = age, y = facial_length, color = dominance_rank, shape = dominance_rank)) +
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.7) +  # Adds a single regression line for all data
  labs(
    title = "Facial Height (mm)",
    x = "Age (years)",
    y = "Facial Height (mm)",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensure the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  # Adjusts the legend position to the bottom
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and pale green for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a5

#Ratio analyses####
# Prepare data for facial width to height ratio
df_ratio_facial <- df_filtered %>%
  filter(!is.na(facial_width) & !is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(facial_ratio = facial_width / facial_length)

df_ratio_facial <- df_ratio_facial %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_ratio))

#Count sample size for this analysis
table(df_ratio_facial$individual, df_ratio_facial$dominance_rank)

#Check assumptions
hist(df_ratio_facial$facial_ratio)
shapiro.test(df_ratio_facial$facial_ratio) #p=0.90

#Fit mixed-effect model
facial_ratio_model <- lme4::lmer(facial_ratio ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_ratio_facial)

#Diagnostics and model output
summary(facial_ratio_model)
drop1(facial_ratio_model, test= "Chisq")
r.squaredGLMM(facial_ratio_model)

# Residuals for facial ratio model
res_facial_ratio_model <- simulateResiduals(facial_ratio_model)
plot(res_facial_ratio_model)
testUniformity(res_facial_ratio_model)
testDispersion(res_facial_ratio_model)
simulationOutput <- simulateResiduals(fittedModel = facial_ratio_model)

# Performing Kolmogorov-Smirnov test for normality
ks_test_ratio <- testUniformity(simulationOutput)
print(ks_test_ratio) 

# Checking for overdispersion
dispersion_test_ratio <- testDispersion(simulationOutput)
print(dispersion_test_ratio) 

#test outliers
testOutliers(res_facial_ratio_model)

# Plot for facial width to height ratio
p6 <- ggplot(df_ratio_facial, aes(x = dominance_rank, y = facial_ratio, fill = dominance_rank)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +  # Use geom_jitter
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.3) +
  labs(y = NULL, x = NULL, title = "Facial Width-to-Height Ratio") +  # Remove y-axis title
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  ) +
  scale_fill_manual(values = c("Alpha" = "#6A0DAD", "Subordinate" = "#A9DFBF")) +
  guides(fill = FALSE)  # Remove legend
p6


# Scatter plot showing facial length vs age by dominance rank
a6 <- ggplot(df_ratio_facial, aes(x = age, y = facial_ratio, color = dominance_rank, shape = dominance_rank)) +
  geom_point(size = 3, alpha = 0.6) +  # Adds points with color and shape differentiation
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black", alpha = 0.7) +  # Adds a single regression line for all data
  labs(
    title = "Facial Width-to-Height Ratio",
    x = "Age (years)",
    y = "FWHR",
    color = "Dominance Rank",
    shape = "Dominance Rank"  # Ensures the legend for shape is also shown
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centers the title
    legend.position = "bottom"  # Adjusts the legend position to the bottom
  ) +
  scale_color_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF")) +  # Deep indigo and pale green for contrast
  scale_shape_manual(values = c("Alpha" = 17, "Subordinate" = 15))  # Different shapes for each rank
a6

#Create plots for scrotum and facial widths####
library(patchwork)

dominance_plot <- (p2 | p3 | p4) / 
  (p5 | p6 | p1)
dominance_plot

# Saving the plot as a high-resolution PNG for digital use
ggsave("dominance_plot.png", plot = dominance_plot, width = 8, height = 6, dpi = 300)

# Remove legends from individual plots if not already done
a1 <- a1 + theme(legend.position = "none") 
a2 <- a2 + theme(legend.position = "none")
a3 <- a3 + theme(legend.position = "none")
a4 <- a4 + theme(legend.position = "none") 
a5 <- a5 + theme(legend.position = "none") 
a6 <- a6 + theme(legend.position = "none") 

# Combine the plots
plot2 <- (a2 | a3 | a4) /
  (a5 | a6 | a1)

# Add a single legend
final_age_plot <- plot2 + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
print(final_age_plot)

# Saving the plot as a high-resolution PNG for digital use
ggsave("Age_plot.png", plot = final_age_plot, width = 8, height = 6, dpi = 300)

#Facial analyses With outlier (LE). Used for supplemental material####
#Facial width analysis
# Data Cleaning
df_facial_widthL <- df %>%
  filter(!is.na(facial_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))
df_facial_widthL <- df_facial_widthL %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_width))

#Check Assumptions with Diagnostics
hist(df_facial_widthL$facial_width) #looks v normal
shapiro.test(df_facial_widthL$facial_width) # normal p=0.85

# count sample size for this analysis
table(df_facial_widthL$individual, df_facial_widthL$dominance_rank)

# Fit mixed-effects model for facial length
facial_width_modell <- lme4::lmer(facial_width ~ dominance_rank + body_length + age + (1 | individual), data = df_facial_widthL)

#Diagnostics and model output
summary(facial_width_modell)
drop1(facial_width_modell, test="Chisq")
r.squaredGLMM(facial_width_modell)

# Simulate residuals and plot diagnostics
res_facial_width_model <- simulateResiduals(facial_width_model)
plot(res_facial_width_model) 
plot(allEffects(facial_width_model)) 

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_fw <- testUniformity(res_facial_width_model)
print(ks_test_result_fw) 

# Perform dispersion test
dispersion_test_result_fw <- testDispersion(res_facial_width_model)
print(dispersion_test_result_fw) 

#test outliers
testOutliers(res_facial_width_model)

#Muzzle width analysis
#Data Cleaning
df_muzzlel <- df %>%
  filter(!is.na(muzzle_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_muzzlel <- df_muzzlel %>%
  filter(complete.cases(dominance_rank, body_length, age, muzzle_width))

#Check assumptions
hist(df_muzzlel$muzzle_width) #looks normal
shapiro.test(df_muzzlel$muzzle_width) #normal p=0.07

#Fix mixed-effects model
muzzle_width_modell <- lme4::lmer(muzzle_width ~ dominance_rank + body_length + age + (1 | individual), data = df_muzzlel)

#Diagnostics and model output
summary(muzzle_width_modell)
drop1(muzzle_width_modell, test="Chisq")
r.squaredGLMM(muzzle_width_modell)

#Model Diagnostics
res_muzzle_width_model <- simulateResiduals(muzzle_width_modell)
plot(allEffects(muzzle_width_modell))

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_muzzle <- testUniformity(res_muzzle_width_modell)
print(ks_test_result_muzzle)

# Perform dispersion test
dispersion_test_result_muzzle <- testDispersion(res_muzzle_width_model)
print(dispersion_test_result_muzzle)

#test for outliers
testOutliers(res_muzzle_width_model)

#Brow width analysis
# Filter out missing values for brow_width
df_brow <- df %>%
  filter(!is.na(brow_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_brow <- df_brow %>%
  filter(complete.cases(dominance_rank, body_length, age, brow_width))

#Check assumptions 
hist(df_brow$brow_width)
shapiro.test(df_brow$brow_width)#not normal.. 

leveneTest(brow_width ~ dominance_rank, data = df_brow) 

# count sample size for this analysis
table(df_brow$individual, df_brow$dominance_rank)

#log transformation for brow
df_brow <- df_brow%>%
  mutate(log_brow_width = log(brow_width))

# Perform Shapiro-Wilk test on the log-transformed data
shapiro_test_log <- shapiro.test(df_brow$log_brow_width)
print(shapiro_test_log) 

#Fit model for log brow width
logbrow_model <- lme4::lmer(log_brow_width ~ dominance_rank + body_length + age + (1 | individual), data = df_brow)

#Diagnostics and model output
summary(logbrow_model)
drop1(logbrow_model, test="Chisq")
r.squaredGLMM(logbrow_model)

#Residuals diagnostics
residuals_sim <- simulateResiduals(fittedModel = logbrow_model)
plot(residuals_sim)

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_brow <- testUniformity(res_brow_width_model)
print(ks_test_result_brow) 

# Perform dispersion test
dispersion_test_result_brow <- testDispersion(res_brow_width_model)
print(dispersion_test_result_brow)

#test for outliers
testOutliers(res_brow_width_model)

#Facial length analysis
# Remove rows with missing values for the columns used in the model
df_lengthl <- df %>%
  filter(!is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual))

df_lengthl <- df_lengthl %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_length))

# Check assumptions
hist(df_lengthl$facial_length)
shapiro.test(df_lengthl$facial_length)  # normal (p=.76)

# count sample size for this analysis
table(df_lengthl$individual, df_lengthl$dominance_rank)

# Fit mixed-effects model for facial length
facial_length_modell <- lme4::lmer(facial_length ~ dominance_rank + body_length + age + (1 | individual), data = df_lengthl)

#Diagnostics and model output
summary(facial_length_modell)
drop1(facial_length_modell, test="Chisq")
r.squaredGLMM(facial_length_modell)

# Simulate residuals and plot diagnostics
res_facial_length_model <- simulateResiduals(facial_length_model)
plot(res_facial_length_model)  # Residual diagnostics
plot(allEffects(facial_length_model))

testOutliers(res_facial_length_model)

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_length <- testUniformity(res_facial_length_model)
print(ks_test_result_length) 

# Perform dispersion test
dispersion_test_result_length <- testDispersion(res_facial_length_model)
print(dispersion_test_result_length)

#facial width-to-height ratio analyses
# Prepare data for facial width to height ratio
df_ratio_faciall <- df %>%
  filter(!is.na(facial_width) & !is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(facial_ratio = facial_width / facial_length)

df_ratio_faciall <- df_ratio_faciall %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_ratio))

# Check assumptions
hist(df_ratio_faciall$facial_ratio)
shapiro.test(df_ratio_faciall$facial_ratio) #p=0.98

# count sample size for this analysis
table(df_ratio_faciall$individual, df_ratio_faciall$dominance_rank)

#Fit model for facial ratio
facial_ratio_modell <- lme4::lmer(facial_ratio ~ dominance_rank + body_length + age + (1 | individual), data = df_ratio_faciall)

#Diagnotics and model output
summary(facial_ratio_modell)
drop1(facial_ratio_modell, test= "Chisq")
r.squaredGLMM(facial_ratio_modell)

# Residuals for facial ratio model
res_facial_ratio_model <- simulateResiduals(facial_ratio_model)
plot(res_facial_ratio_model)
testUniformity(res_facial_ratio_model)
testDispersion(res_facial_ratio_model)
simulationOutput <- simulateResiduals(fittedModel = facial_ratio_model)

# Performing Kolmogorov-Smirnov test 
ks_test_ratio <- testUniformity(simulationOutput)
print(ks_test_ratio)

# Checking for overdispersion
dispersion_test_ratio <- testDispersion(simulationOutput)
print(dispersion_test_ratio)

#IOC Calculations####
# Convert all columns to character first
iior <- iior %>% mutate(across(everything(), as.character))

# Replace '#VALUE!' with NA
iior_clean <- iior %>% mutate(across(everything(), ~na_if(.x, "#VALUE!")))

# View cleaned data
print(iior_clean)
head(iior_clean)
iior_clean <- iior_clean %>% select(1:8)
str(iior_clean)

# Convert columns 3 to 8 to numeric
iior_clean <- iior_clean %>% mutate(across(3:8, as.numeric))

iior_long <- iior_clean %>%
  pivot_longer(cols = brow:body, 
               names_to = "measurement", 
               values_to = "value")

#reorganize
ior <-iior_long
head(ior)

# Brow Width IOC 
# Filter the dataset for the 'brow' measurement
brow_data <- ior %>% filter(measurement == "brow")

# Pivot the data to a wider format, ensuring no duplicate warnings
brow_data_wide <- brow_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
brow_data_wide <- brow_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

# Check the cleaned data
print(brow_data_wide)

brow_data_wide <- brow_data_wide %>%
  slice_head(n = 50)

# Convert to matrix for ICC calculation
brow_matrix <- as.matrix(brow_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_brow <- icc(brow_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_brow)

# Muzzle Width IOC
# Filter the dataset for the 'muzzle' measurement
muzzle_data <- ior %>% filter(measurement == "muzzle")

# Pivot the data to a wider format, ensuring no duplicate warnings
muzzle_data_wide <- muzzle_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
muzzle_data_wide <- muzzle_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

muzzle_data_wide <- muzzle_data_wide %>%
  slice_head(n = 50)

# Convert to matrix for ICC calculation
muzzle_matrix <- as.matrix(muzzle_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_muzzle <- icc(muzzle_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_muzzle)

# Facial Width IOC
# Filter the dataset for the 'facial' measurement
facial_data <- ior %>% filter(measurement == "facial")

# Pivot the data to a wider format, ensuring no duplicate warnings
facial_data_wide <- facial_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
facial_data_wide <- facial_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

facial_data_wide <- facial_data_wide %>%
  slice_head(n = 50)

# Convert to matrix for ICC calculation
facial_matrix <- as.matrix(facial_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_facial <- icc(facial_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_facial)

# Facial Height IOC 
# Filter the dataset for the 'height' measurement
height_data <- ior %>% filter(measurement == "height")

# Pivot the data to a wider format, ensuring no duplicate warnings
height_data_wide <- height_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
height_data_wide <- height_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

height_data_wide <- height_data_wide %>%
  slice_head(n = 50)

# Convert to matrix for ICC calculation
height_matrix <- as.matrix(height_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_height <- icc(height_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_height)

# Scrotum width IOC
# Filter the dataset for the 'scrotum' measurement
scrotum_data <- ior %>% filter(measurement == "scrotum")

# Pivot the data to a wider format, ensuring no duplicate warnings
scrotum_data_wide <- scrotum_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
scrotum_data_wide <- scrotum_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

scrotum_data_wide <- scrotum_data_wide %>%
  slice_head(n = 50)

# Convert to matrix for ICC calculation
scrotum_matrix <- as.matrix(scrotum_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_scrotum <- icc(scrotum_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_scrotum)

# Body length IOC
# Filter the dataset for the 'body' measurement
body_data <- ior %>% filter(measurement == "body")

# Pivot the data to a wider format, ensuring no duplicate warnings
body_data_wide <- body_data %>%
  pivot_wider(names_from = id, values_from = value, values_fn = list)

# Remove list-columns by unnesting and handling rows with complete observations only
body_data_wide <- body_data_wide %>%
  unnest(c(nc, MM, eh)) %>%
  filter(!is.na(nc) & !is.na(MM) & !is.na(eh)) %>%
  filter(nc != 0 & MM != 0 & eh != 0)

body_data_wide <- body_data_wide %>%
  slice_head(n = 50)


# Convert to matrix for ICC calculation
body_matrix <- as.matrix(body_data_wide[, c("nc", "MM", "eh")])

# Calculate ICC
icc_body <- icc(body_matrix, model = "twoway", type = "agreement", unit = "single")

# Print the results
print(icc_body)


# Convert dominance rank to a binary factor for logistic regression
df$dominance_rank_binary <- ifelse(df$dominance_rank == "Alpha", 1, 0)

# Fit logistic regression model
dominance_model <- glm(dominance_rank_binary ~ body_length, data = df, family = "binomial")

# Check the summary of the model to see if body length is a significant predictor
summary(dominance_model)


#Analyses with Laser type included in models####

scrotum_model1 <- lme4::lmer(scrotum_width ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_scrotum_width)
# Diagnostics and model output
summary(scrotum_model1)
drop1(scrotum_model1, test ="Chisq")
r.squaredGLMM(scrotum_model1)

# Simulate and plot residuals
res_scrotum_model1 <- simulateResiduals(scrotum_model1)
plot(res_scrotum_model1)  # Residual diagnostics #LOOKS GOOD
plot(allEffects(scrotum_model1))

# Fit mixed-effects model for facial length
facial_width_model1 <- lme4::lmer(facial_width ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_facial_width)

#Diagnostics and model output
summary(facial_width_model1)
drop1(facial_width_model1, test ="Chisq")
r.squaredGLMM(facial_width_model1)

# Simulate residuals and plot diagnostics
res_facial_width_model1 <- simulateResiduals(facial_width_model1)
plot(res_facial_width_model)  # looks fine
plot(allEffects(facial_width_model)) #looks fine


#Fix mixed-effects model
muzzle_width_model1 <- lme4::lmer(muzzle_width ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_muzzle)

#Diagnostics and model output
summary(muzzle_width_model1)
drop1(muzzle_width_model1, test ="Chisq")
r.squaredGLMM(muzzle_width_model1)

#Model Diagnostics
res_muzzle_width_model1 <- simulateResiduals(muzzle_width_model1)
plot(res_muzzle_width_model1)  
plot(allEffects(muzzle_width_model1))

#Fit mixed-effects model
brow_model1 <- lme4::lmer(brow_width ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_brow)

#Diagnostics and model output
summary(brow_model1)
drop1(brow_model1, test="Chisq")
r.squaredGLMM(brow_model1)

#Residuals diagnostics
residuals_sim1 <- simulateResiduals(fittedModel = brow_model1)
plot(residuals_sim1)
plot(allEffects(brow_model1))

# Fit mixed-effects model for facial length
facial_length_model1 <- lme4::lmer(facial_length ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_length)

#Diagnostics and model output
summary(facial_length_model1)
drop1(facial_length_model1, test= "Chisq")
r.squaredGLMM(facial_length_model1)

# Simulate residuals and plot diagnostics
res_facial_length_model1 <- simulateResiduals(facial_length_model1)
plot(res_facial_length_model1)  # Residual diagnostics
plot(allEffects(facial_length_model1))

#Fit mixed-effect model
facial_ratio_model1 <- lme4::lmer(facial_ratio ~ dominance_rank + body_length + age + laser_type + (1 | individual), data = df_ratio_facial)

#Diagnostics and model output
summary(facial_ratio_model1)
drop1(facial_ratio_model1, test= "Chisq")
r.squaredGLMM(facial_ratio_model1)

# Residuals for facial ratio model
res_facial_ratio_model1 <- simulateResiduals(facial_ratio_model1)
plot(res_facial_ratio_model1)
testUniformity(res_facial_ratio_model1)
testDispersion(res_facial_ratio_model1)

