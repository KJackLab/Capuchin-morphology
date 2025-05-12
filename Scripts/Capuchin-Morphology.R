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
library(emmeans)
library(ggeffects)

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

# Calculate mean scrotum width per individual
df_scrotum_summary <- df_scrotum_width %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_scrotum_width = mean(scrotum_width, na.rm = TRUE),
    sd_scrotum_width = sd(scrotum_width, na.rm = TRUE),
    n = n(),
    se_scrotum_width = sd_scrotum_width / sqrt(n)
  ) %>%
  ungroup()

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
res_scrotum_model <- simulateResiduals(scrotum_model)
plot(res_scrotum_model)  # Residual diagnostics #LOOKS GOOD
plot(allEffects(scrotum_model))

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

# Get model predictions (estimated marginal means)
emm_scrotum <- emmeans(scrotum_model, ~ dominance_rank)

# Convert to dataframe for plotting
emm_scrotum_df <- as.data.frame(emm_scrotum)

set.seed(123)

df_scrotum_summary <- df_scrotum_summary %>%
  mutate(
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")),
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)  # slight jitter
  )

p_scrotum <- ggplot() +
  # Points for individual monkeys
  geom_point(data = df_scrotum_summary,
             aes(x = x_pos, y = mean_scrotum_width),
             size = 3, alpha = 0.7, shape = 21, fill = "white", color = "black") +
  
  # Error bars for individuals
  geom_errorbar(data = df_scrotum_summary,
                aes(x = x_pos, ymin = mean_scrotum_width - se_scrotum_width, ymax = mean_scrotum_width + se_scrotum_width),
                width = 0.05, alpha = 0.5) +
  
  # Model prediction points (no jitter)
  geom_point(data = emm_scrotum_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_scrotum_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +  
  
  labs(x = NULL, y = NULL, title = "Scrotum Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

p_scrotum

#data for figure 5- age plot
# Create plotting dataset with LE included
df_scrotum_age_plot <- df %>%
  filter(!is.na(scrotum_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_scrotum_width = mean(scrotum_width, na.rm = TRUE),
    sd_scrotum_width = sd(scrotum_width, na.rm = TRUE),
    se_scrotum_width = sd_scrotum_width / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")  # Mark LE separately
  )


scrotum_age_pred <- ggpredict(scrotum_model, terms = "age [all]")


p_scrotum_age <- ggplot() +
  geom_ribbon(data = scrotum_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = scrotum_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  geom_point(data = df_scrotum_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_scrotum_width, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +  # Black outline, solid fill
  
  # LE point separately (red fill)
  geom_point(data = df_scrotum_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_scrotum_width),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  geom_errorbar(data = df_scrotum_age_plot,
                aes(x = mean_age, ymin = mean_scrotum_width - se_scrotum_width, ymax = mean_scrotum_width + se_scrotum_width),
                width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  labs(x = "Age (years)", y = NULL, title = "Scrotum Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # Remove legend
  )

p_scrotum_age

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

lme4::ranef(facial_width_model)
ranefs <- ranef(facial_width_model)$individual
ranefs$ID <- rownames(ranefs)
ggplot(ranefs, aes(x = reorder(ID, `(Intercept)`), y = `(Intercept)`)) +
  geom_point() +
  coord_flip() +
  labs(y = "Random Intercept (BLUP)", x = "Individual", title = "Individual Deviations in Facial Width")


# Perform dispersion test
dispersion_test_result_fw <- testDispersion(res_facial_width_model)
# Print dispersion test result
print(dispersion_test_result_fw) #fine p=0.63

#test outliers
testOutliers(res_facial_width_model)

#model predictions for plot (without outlier)
emm_facial_width <- emmeans(facial_width_model, ~ dominance_rank)
emm_facial_width_df <- as.data.frame(emm_facial_width)

#data for plot (with outlier)
df_facial_width_plot <- dff %>%
  filter(!is.na(facial_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual)) %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_width)) %>%
  mutate(
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")),
    facial_width = as.numeric(facial_width)  
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_width = mean(facial_width, na.rm = TRUE),
    sd_facial_width = sd(facial_width, na.rm = TRUE),
    n = n(),
    se_facial_width = sd_facial_width / sqrt(n)
  ) %>%
  ungroup()

#Add LE points as red
df_facial_width_plot <- df_facial_width_plot %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )


#plot preparation
set.seed(123)

df_facial_width_plot <- df_facial_width_plot %>%
  mutate(
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)
  )


#plot of facial width by dominance
p_facial_width <- ggplot() +
  # Individual points (color by LE vs Other)
  geom_point(data = df_facial_width_plot,
             aes(x = x_pos, y = mean_facial_width, fill = is_LE),
             size = 3, alpha = 0.7, shape = 21, color = "black") +  # Shape 21 = hollow circle
  geom_errorbar(data = df_facial_width_plot,
                aes(x = x_pos, ymin = mean_facial_width - se_facial_width, ymax = mean_facial_width + se_facial_width),
                width = 0.05, alpha = 0.5) +
  
  # Model predictions
  geom_point(data = emm_facial_width_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_facial_width_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_fill_manual(values = c("Other" = "white", "LE" = "red")) + 
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +
  
  labs(x = NULL, y = NULL, title = "Facial Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")  

p_facial_width

#data for figure 5-age plot
# Create plotting dataset with LE included
# Create plotting dataset with LE included
df_facial_width_age_plot <- dff %>%
  filter(!is.na(facial_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(
    facial_width = as.numeric(facial_width),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")) 
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_width = mean(facial_width, na.rm = TRUE),
    sd_facial_width = sd(facial_width, na.rm = TRUE),
    se_facial_width = sd_facial_width / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )


#to plot model predictions
facial_width_age_pred <- ggpredict(facial_width_model, terms = "age [all]")

#plot for age vs facial width
p_facial_width_age <- ggplot() +
  # Model prediction line + CI shading
  geom_ribbon(data = facial_width_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = facial_width_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  
  # Points for non-LE individuals, correctly mapped FILL
  geom_point(data = df_facial_width_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_facial_width, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +
  
  # LE point separately (solid red)
  geom_point(data = df_facial_width_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_facial_width),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  
  # Error bars for individuals
  geom_errorbar(data = df_facial_width_age_plot,
                aes(x = mean_age, ymin = mean_facial_width - se_facial_width, ymax = mean_facial_width + se_facial_width),
                width = 0.2, alpha = 0.5) +
  
  # Custom fill scale
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  
  labs(x = "Age (years)", y = NULL, title = "Facial Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # No legend
  )

p_facial_width_age

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

#model predictions for plot
emm_muzzle_width <- emmeans(muzzle_width_model, ~ dominance_rank)
emm_muzzle_width_df <- as.data.frame(emm_muzzle_width)

# Data for plotting (with LE)
df_muzzle_width_plot <- dff %>%
  filter(!is.na(muzzle_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual)) %>%
  filter(complete.cases(dominance_rank, body_length, age, muzzle_width)) %>%
  mutate(
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")),
    muzzle_width = as.numeric(muzzle_width)  
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_muzzle_width = mean(muzzle_width, na.rm = TRUE),
    sd_muzzle_width = sd(muzzle_width, na.rm = TRUE),
    n = n(),
    se_muzzle_width = sd_muzzle_width / sqrt(n)
  ) %>%
  ungroup()

# Add LE points in red
df_muzzle_width_plot <- df_muzzle_width_plot %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

# Prepare jittered x-positions
set.seed(123)

df_muzzle_width_plot <- df_muzzle_width_plot %>%
  mutate(
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)
  )

#muzzle width by dominance rank
p_muzzle_width <- ggplot() +
  # Individual points (color by LE vs Other)
  geom_point(data = df_muzzle_width_plot,
             aes(x = x_pos, y = mean_muzzle_width, fill = is_LE),
             size = 3, alpha = 0.7, shape = 21, color = "black") +
  
  # Individual error bars
  geom_errorbar(data = df_muzzle_width_plot,
                aes(x = x_pos, ymin = mean_muzzle_width - se_muzzle_width, ymax = mean_muzzle_width + se_muzzle_width),
                width = 0.05, alpha = 0.5) +
  
  # Model predictions
  geom_point(data = emm_muzzle_width_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_muzzle_width_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_fill_manual(values = c("Other" = "white", "LE" = "red")) +
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +
  
  labs(x = NULL, y = NULL, title = "Muzzle Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # Optional: remove legend

p_muzzle_width

#plot for figure 5 age vs muzzle
# Create plotting dataset with LE included
df_muzzle_width_age_plot <- dff %>%
  filter(!is.na(muzzle_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(
    muzzle_width = as.numeric(muzzle_width),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate"))
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_muzzle_width = mean(muzzle_width, na.rm = TRUE),
    sd_muzzle_width = sd(muzzle_width, na.rm = TRUE),
    se_muzzle_width = sd_muzzle_width / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )


muzzle_width_age_pred <- ggpredict(muzzle_width_model, terms = "age [all]")


p_muzzle_width_age <- ggplot() +
  # Model prediction line + CI shading
  geom_ribbon(data = muzzle_width_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = muzzle_width_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  
  # Points for non-LE individuals
  geom_point(data = df_muzzle_width_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_muzzle_width, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +
  
  # LE point separately (solid red)
  geom_point(data = df_muzzle_width_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_muzzle_width),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  
  # Error bars for individuals
  geom_errorbar(data = df_muzzle_width_age_plot,
                aes(x = mean_age, ymin = mean_muzzle_width - se_muzzle_width, ymax = mean_muzzle_width + se_muzzle_width),
                width = 0.2, alpha = 0.5) +
  
  # Custom fill scale
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  
  labs(x = "Age (years)", y = NULL, title = "Muzzle Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

p_muzzle_width_age


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
res_brow_width_model<- simulateResiduals(fittedModel = brow_model)
plot(res_brow_width_model)

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

#model estimates for plot
emm_brow_width <- emmeans(brow_model, ~ dominance_rank)
emm_brow_width_df <- as.data.frame(emm_brow_width)

#data for plot
df_brow_width_plot <- dff %>%
  filter(!is.na(brow_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual)) %>%
  filter(complete.cases(dominance_rank, body_length, age, brow_width)) %>%
  mutate(
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")),
    brow_width = as.numeric(brow_width)  
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_brow_width = mean(brow_width, na.rm = TRUE),
    sd_brow_width = sd(brow_width, na.rm = TRUE),
    n = n(),
    se_brow_width = sd_brow_width / sqrt(n)
  ) %>%
  ungroup()

# Add LE points
df_brow_width_plot <- df_brow_width_plot %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

# Jitter x-axis
set.seed(123)

df_brow_width_plot <- df_brow_width_plot %>%
  mutate(
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)
  )

#brow width vs dominance plot
p_brow_width <- ggplot() +
  # Individual points (color by LE)
  geom_point(data = df_brow_width_plot,
             aes(x = x_pos, y = mean_brow_width, fill = is_LE),
             size = 3, alpha = 0.7, shape = 21, color = "black") +
  
  # Individual error bars
  geom_errorbar(data = df_brow_width_plot,
                aes(x = x_pos, ymin = mean_brow_width - se_brow_width, ymax = mean_brow_width + se_brow_width),
                width = 0.05, alpha = 0.5) +
  
  # Model predictions
  geom_point(data = emm_brow_width_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_brow_width_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_fill_manual(values = c("Other" = "white", "LE" = "red")) +
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +
  
  labs(x = NULL, y = NULL, title = "Brow Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p_brow_width

#plot for figure 5 age vs muzzle
# Create plotting dataset with LE included
df_brow_width_age_plot <- dff %>%
  filter(!is.na(brow_width) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(
    brow_width = as.numeric(brow_width),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate"))
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_brow_width = mean(brow_width, na.rm = TRUE),
    sd_brow_width = sd(brow_width, na.rm = TRUE),
    se_brow_width = sd_brow_width / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

#plot model predictions
brow_width_age_pred <- ggpredict(brow_model, terms = "age [all]")

p_brow_width_age <- ggplot() +
  # Model prediction line + CI shading
  geom_ribbon(data = brow_width_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = brow_width_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  
  # Points for non-LE individuals
  geom_point(data = df_brow_width_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_brow_width, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +
  
  # LE point separately (solid red)
  geom_point(data = df_brow_width_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_brow_width),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  
  # Error bars for individuals
  geom_errorbar(data = df_brow_width_age_plot,
                aes(x = mean_age, ymin = mean_brow_width - se_brow_width, ymax = mean_brow_width + se_brow_width),
                width = 0.2, alpha = 0.5) +
  
  # Custom fill scale
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  
  labs(x = "Age (years)", y = NULL, title = "Brow Width (mm)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

p_brow_width_age



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

#model predictions for plot
emm_facial_length <- emmeans(facial_length_model, ~ dominance_rank)
emm_facial_length_df <- as.data.frame(emm_facial_length)

#data for plot
df_facial_length_plot <- dff %>%
  filter(!is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual)) %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_length)) %>%
  mutate(
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate")),
    facial_length = as.numeric(facial_length)  
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_length = mean(facial_length, na.rm = TRUE),
    sd_facial_length = sd(facial_length, na.rm = TRUE),
    n = n(),
    se_facial_length = sd_facial_length / sqrt(n)
  ) %>%
  ungroup()

# Add LE indicator
df_facial_length_plot <- df_facial_length_plot %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

# Jittered x positions
set.seed(123)

df_facial_length_plot <- df_facial_length_plot %>%
  mutate(
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)
  )

p_facial_length <- ggplot() +
  # First: plot all the OTHER individuals (white dots)
  geom_point(data = df_facial_length_plot %>% filter(is_LE == "Other"),
             aes(x = x_pos, y = mean_facial_length),
             size = 3, alpha = 0.7, shape = 21, fill = "white", color = "black") +
  
  # Then: plot LE separately (red dot)
  geom_point(data = df_facial_length_plot %>% filter(is_LE == "LE"),
             aes(x = x_pos, y = mean_facial_length),
             size = 3, alpha = 0.9, shape = 21, fill = "red", color = "black") +
  
  # Individual error bars (for everyone)
  geom_errorbar(data = df_facial_length_plot,
                aes(x = x_pos, ymin = mean_facial_length - se_facial_length, ymax = mean_facial_length + se_facial_length),
                width = 0.05, alpha = 0.5) +
  
  # Model predictions
  geom_point(data = emm_facial_length_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_facial_length_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +
  
  labs(x = NULL, y = NULL, title = "Facial Height (mm)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p_facial_length

#plot for figure 5- age vs facial height
# Create plotting dataset with LE included
df_facial_length_age_plot <- dff %>%
  filter(!is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(
    facial_length = as.numeric(facial_length),
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate"))
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_length = mean(facial_length, na.rm = TRUE),
    sd_facial_length = sd(facial_length, na.rm = TRUE),
    se_facial_length = sd_facial_length / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

#plot model predictions
facial_length_age_pred <- ggpredict(facial_length_model, terms = "age [all]")

p_facial_length_age <- ggplot() +
  # Model prediction line + CI shading
  geom_ribbon(data = facial_length_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = facial_length_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  
  # Points for non-LE individuals
  geom_point(data = df_facial_length_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_facial_length, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +
  
  # LE point separately (solid red)
  geom_point(data = df_facial_length_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_facial_length),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  
  # Error bars for individuals
  geom_errorbar(data = df_facial_length_age_plot,
                aes(x = mean_age, ymin = mean_facial_length - se_facial_length, ymax = mean_facial_length + se_facial_length),
                width = 0.2, alpha = 0.5) +
  
  # Custom fill scale
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  
  labs(x = "Age (years)", y = NULL, title = "Facial Height (mm)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

p_facial_length_age

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


# model predictions for plot
emm_facial_ratio <- emmeans(facial_ratio_model, ~ dominance_rank)
emm_facial_ratio_df <- as.data.frame(emm_facial_ratio)

# Data for facial width-to-height ratio (fWHR) plot
df_facial_ratio_plot <- dff %>%
  mutate(
    facial_ratio = as.numeric(facial_width) / as.numeric(facial_length)  # Define facial ratio first
  ) %>%
  filter(!is.na(facial_ratio) & !is.na(dominance_rank) & !is.na(age) & !is.na(Temp) & !is.na(individual)) %>%
  filter(complete.cases(dominance_rank, body_length, age, facial_width, facial_length)) %>%
  mutate(
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate"))
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_ratio = mean(facial_ratio, na.rm = TRUE),
    sd_facial_ratio = sd(facial_ratio, na.rm = TRUE),
    n = n(),
    se_facial_ratio = sd_facial_ratio / sqrt(n)
  ) %>%
  ungroup()

# Add LE indicator
df_facial_ratio_plot <- df_facial_ratio_plot %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

# Jitter x positions
set.seed(123)

df_facial_ratio_plot <- df_facial_ratio_plot %>%
  mutate(
    x_pos = ifelse(dominance_rank == "Alpha", 1, 1.5) + runif(n(), min = -0.1, max = 0.1)
  )

p_facial_ratio <- ggplot() +
  # First: plot all the OTHER individuals
  geom_point(data = df_facial_ratio_plot %>% filter(is_LE == "Other"),
             aes(x = x_pos, y = mean_facial_ratio),
             size = 3, alpha = 0.7, shape = 21, fill = "white", color = "black") +
  
  # Then: plot LE separately
  geom_point(data = df_facial_ratio_plot %>% filter(is_LE == "LE"),
             aes(x = x_pos, y = mean_facial_ratio),
             size = 3, alpha = 0.9, shape = 21, fill = "red", color = "black") +
  
  # Individual error bars
  geom_errorbar(data = df_facial_ratio_plot,
                aes(x = x_pos, ymin = mean_facial_ratio - se_facial_ratio, ymax = mean_facial_ratio + se_facial_ratio),
                width = 0.05, alpha = 0.5) +
  
  # Model predictions
  geom_point(data = emm_facial_ratio_df,
             aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), y = emmean),
             size = 5, shape = 23, fill = "black", color = "black") +
  
  geom_errorbar(data = emm_facial_ratio_df,
                aes(x = ifelse(dominance_rank == "Alpha", 1, 1.5), ymin = lower.CL, ymax = upper.CL),
                width = 0.2, color = "black") +
  
  scale_x_continuous(breaks = c(1, 1.5),
                     labels = c("Alpha", "Subordinate"),
                     limits = c(0.8, 1.7)) +
  
  labs(x = NULL, y = NULL, title = "Facial Width-to-Height Ratio (fWHR)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size =15),
        legend.position = "none")

p_facial_ratio

#plot for figure 5- age vs rati
# Create plotting dataset with LE included and define ratio
df_facial_ratio_age_plot <- dff %>%
  filter(!is.na(facial_width) & !is.na(facial_length) & !is.na(dominance_rank) & !is.na(age) & !is.na(individual)) %>%
  mutate(
    facial_width = as.numeric(facial_width),
    facial_length = as.numeric(facial_length),
    facial_ratio = facial_width / facial_length,  # Define facial ratio here
    dominance_rank = case_when(
      dominance_rank == "a" ~ "Alpha",
      dominance_rank == "s" ~ "Subordinate",
      TRUE ~ dominance_rank
    ),
    dominance_rank = factor(dominance_rank, levels = c("Alpha", "Subordinate"))
  ) %>%
  group_by(individual, dominance_rank) %>%
  summarise(
    mean_facial_ratio = mean(facial_ratio, na.rm = TRUE),
    sd_facial_ratio = sd(facial_ratio, na.rm = TRUE),
    se_facial_ratio = sd_facial_ratio / sqrt(n()),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    is_LE = ifelse(individual == "LE", "LE", "Other")
  )

facial_ratio_age_pred <- ggpredict(facial_ratio_model, terms = "age [all]")

p_facial_ratio_age <- ggplot() +
  # Model prediction line + CI shading
  geom_ribbon(data = facial_ratio_age_pred, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "grey") +
  geom_line(data = facial_ratio_age_pred, aes(x = x, y = predicted), size = 1, color = "black") +
  
  # Points for non-LE individuals
  geom_point(data = df_facial_ratio_age_plot %>% filter(is_LE == "Other"),
             aes(x = mean_age, y = mean_facial_ratio, fill = dominance_rank),
             size = 3, shape = 21, color = "black", alpha = 0.9) +
  
  # LE point separately (solid red)
  geom_point(data = df_facial_ratio_age_plot %>% filter(is_LE == "LE"),
             aes(x = mean_age, y = mean_facial_ratio),
             size = 3, shape = 21, fill = "red", color = "black", alpha = 0.9) +
  
  # Error bars for individuals
  geom_errorbar(data = df_facial_ratio_age_plot,
                aes(x = mean_age, ymin = mean_facial_ratio - se_facial_ratio, ymax = mean_facial_ratio + se_facial_ratio),
                width = 0.2, alpha = 0.5) +
  
  # Custom fill scale
  scale_fill_manual(values = c("Alpha" = "#4B0082", "Subordinate" = "#A9DFBF"), name = "Dominance Rank") +
  
  labs(x = "Age (years)", y = NULL, title = "Facial Width-to-Height Ratio (fWHR)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none"
  )

p_facial_ratio_age


#Create plots for scrotum and facial widths####
library(patchwork)

dominance_plot <- (p_facial_width | p_muzzle_width | p_brow_width) / 
  (p_facial_length | p_facial_ratio | p_scrotum)

dominance_plot

# Saving the plot as a high-resolution PNG for digital use
ggsave("~/Desktop/dominance_plot.png", plot = dominance_plot, width = 10, height = 6, dpi = 300)

# Remove legends individually if needed (but you already set legend.position = "none" in each plot earlier)
p_scrotum_age <- p_scrotum_age + theme(legend.position = "none")
p_facial_width_age <- p_facial_width_age + theme(legend.position = "none")
p_muzzle_width_age <- p_muzzle_width_age + theme(legend.position = "none")
p_brow_width_age <- p_brow_width_age + theme(legend.position = "none")
p_facial_length_age <- p_facial_length_age + theme(legend.position = "none")
p_facial_ratio_age <- p_facial_ratio_age + theme(legend.position = "none")

# Combine the plots
plot2 <- (p_facial_width_age | p_muzzle_width_age | p_brow_width_age) /
  (p_facial_length_age | p_facial_ratio_age | p_scrotum_age)

# Add a single combined guide (legend), and put it at the bottom
final_age_plot <- plot2 + plot_layout(guides = 'collect') & theme(legend.position = "bottom")

# Print the final figure
print(final_age_plot)

# Saving the plot as a high-resolution PNG for digital use
ggsave("~/Desktop/Age_plot.png", plot = final_age_plot, width = 10, height = 6, dpi = 300)


getwd()




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
res_facial_width_modell <- simulateResiduals(facial_width_modell)
plot(res_facial_width_modell) 
plot(allEffects(facial_width_modell)) 

# Perform Kolmogorov-Smirnov (KS) test
ks_test_result_fw <- testUniformity(res_facial_width_modell)
print(ks_test_result_fw) 

# Perform dispersion test
dispersion_test_result_fw <- testDispersion(res_facial_width_modell)
print(dispersion_test_result_fw) 

#test outliers
testOutliers(res_facial_width_modell)

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



