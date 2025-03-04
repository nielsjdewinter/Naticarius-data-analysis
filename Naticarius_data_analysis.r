# Data analysis routine for the manuscript titled "Clumped isotope compositions of mid-Holocene marine gastropod shells"
# Authors of manuscript: Tingting Yua, Sha Lia, Niels J. de Winter b, c, Dongxv Zhanga, Naihua Xuec, d, Ye Xua, Haichun Zhanga, Xiaoqiao Wane, Bo Wanga
# Author of the data analysis code: Niels J. de Winter
# E-mail: n.j.de.winter@vu.nl

# Load necessary libraries
library(ggplot2)
library(tidyverse)

# Load the CSV file
data <- read.csv("Naticarius_all_data.csv")

# Display the first few rows of the data
head(data)

# Rename the columns
colnames(data) <- gsub("^d13C.*", "d13C", colnames(data))
colnames(data) <- gsub("^d18O.*", "d18O", colnames(data))
colnames(data) <- gsub("^D47.*", "D47", colnames(data))

# Filter the data for type "standard"
standard_data <- subset(data, Type == "standard")

# Convert the ID column to Date type
standard_data$Date <- as.POSIXct(standard_data$ID, format = "%Y-%m-%d %H:%M")

# Reshape the data for ggplot
standard_data_long <- standard_data %>%
    pivot_longer(cols = c("d13C", "d18O", "D47"), names_to = "variable", values_to = "value")

# Plot the data using ggplot2 with each standard name having its own panel
ggplot(standard_data_long, aes(x = Date, y = value, color = Easotope.Name)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ Easotope.Name + variable, scales = "free_y", ncol = 3, nrow = 4) +
    labs(x = "Date", y = "Value", title = "Isotope Data over Time for Standard Data") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(aspect.ratio = 1 / 1.5)

# Calculate the standard deviations on all ETH standards
STD_stats <- standard_data_long %>%
    group_by(Easotope.Name, variable) %>%
    summarise(
        mean = mean(value),
        sd = sd(value),
        N = n(),
        se = sd / sqrt(N),
        CL95 = qt(0.975, N - 1) * se
    )

# Calculate the mean standard deviation for D47 based on the standard deviations of the ETH standards
SD_ETH <- sqrt(sum(STD_stats$sd[STD_stats$variable == "D47"] ** 2) / length(unique(STD_stats$Easotope.Name[STD_stats$variable == "D47"])))
data$D47_SD <- SD_ETH # Add ETH standard deviation as column to all data

# ------------------------------------------------------------------------------------------
# Calculate pooled stats for all specimens
# Group data by Sample_name and calculate the mean and the pooled standard deviation for D47
pooled_sample_stats <- data %>%
    subset(Type == "sample") %>%
    group_by(Sample_name) %>%
    summarise(
        N = n(),
        mean_D47 = mean(D47, na.rm = TRUE),
        sd_D47 = sqrt(sum((D47 - mean(D47, na.rm = TRUE)) ** 2, na.rm = TRUE) / (n() - 1)),
        D47_pooled_SD = sqrt(sd_D47 ** 2 + D47_SD ** 2),
        D47_se = D47_pooled_SD / sqrt(N),
        D47_CL95 = qt(0.975, N - 1) * D47_se
    ) %>%
    distinct(Sample_name, .keep_all = TRUE)

# Merge the pooled statistics back into the original data
sample_data <- merge(data, pooled_sample_stats, by = "Sample_name") %>%
    distinct(Sample_name, .keep_all = TRUE)

# Convert specimen mean and 95%CL to temperature
sample_data$T47_sample <- sqrt(10 ** 6 / ((sample_data$D47 - 0.154) / 0.0391)) - 273.15 # Convert D47 to temperature using Anderson equation
sample_data$T47_max_sample <- sqrt(10 ** 6 / ((sample_data$D47 - 2 * sample_data$D47_se - 0.154) / 0.0391)) - 273.15 # Convert D47 - 2 standard errors to temperature using Anderson equation
sample_data$T47_min_sample <- sqrt(10 ** 6 / ((sample_data$D47 + 2 * sample_data$D47_se - 0.154) / 0.0391)) - 273.15 # Convert D47 + 2 standard errors to temperature using Anderson equation
sample_data$T47_95CL_sample <- (sample_data$T47_max_sample - sample_data$T47_min_sample) / 2

# ------------------------------------------------------------------------------------------
# Add calculations of the d18Ow of the water from D47-temperature and d18Oc
# Using the d18Oc-d18Ow-temperature relationship from Grossman and Ku (1986)
# with the VPDB-VSMOW scale correction of ‐0.27 ‰ (Gonfiantini et al., 1995; Dettman et al., 1999)
# Solved for d18Ow following Schöne et al. (2020; see Schmitt et al., 2024)
# Schmitt, K.E., Beuzen-Waller, T., Schmidt, C., Proctor, L., Lindauer, S., Gey, C.J., Pietsch, D., Schöne, B.R., 2024. Melanoides tuberculata and Zootecus insularis gastropod shells provide a snapshot into past hydroclimatic conditions of arid environments: New perspectives from Oman. Palaeogeography, Palaeoclimatology, Palaeoecology 655, 112542. https://doi.org/10.1016/j.palaeo.2024.112542
sample_data$d18Ow <- (19.43 - 4.34 * sample_data$d18O - sample_data$T47_sample) / -4.34
sample_data$d18Ow_min <- (19.43 - 4.34 * sample_data$d18O - sample_data$T47_max_sample) / -4.34
sample_data$d18Ow_max <- (19.43 - 4.34 * sample_data$d18O - sample_data$T47_min_sample) / -4.34
sample_data$d18Ow_95CL <- (sample_data$d18Ow_max - sample_data$d18Ow_min) / 2

# ------------------------------------------------------------------------------------------
# Calculate mean temperature +/- 95%CL for entire dataset
# Calculate the mean and the pooled standard deviation for D47
pooled_data_stats <- data %>%
    subset(Type == "sample") %>%
    summarise(
        N = n(),
        mean_D47 = mean(D47, na.rm = TRUE),
        D47_pooled_SD = sqrt(sum((D47 - mean(D47, na.rm = TRUE)) ** 2, na.rm = TRUE) / (n() - 1)),
        D47_se = D47_pooled_SD / sqrt(N),
        D47_CL95 = qt(0.975, N - 1) * D47_se,
        T47_mean_all = sqrt(10 ** 6 / ((mean_D47 - 0.154) / 0.0391)) - 273.15, # Convert D47 to temperature using Anderson equation
        T47_max = sqrt(10 ** 6 / ((mean_D47 - D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 - 2 standard errors to temperature using Anderson equation
        T47_min = sqrt(10 ** 6 / ((mean_D47 + D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 + 2 standard errors to temperature using Anderson equation
        T47_95CL_all = (T47_max - T47_min) / 2,
        d18Ow_mean_all = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_mean_all) / -4.34, # Calculate d18Ow from mean d18Oc and mean T47
        d18Ow_max = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_max) / -4.34, # Calculate d18Ow from mean d18Oc and T47 - 2 standard errors
        d18Ow_min = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_min) / -4.34, # Calculate d18Ow from mean d18Oc and T47 + 2 standard errors
        d18Ow_95CL_all = (d18Ow_max - d18Ow_min) / 2
    )

# Merge the pooled statistics back into the original data
sample_data$T47_mean_all <- pooled_data_stats$T47_mean_all
sample_data$T47_95CL_all <- pooled_data_stats$T47_95CL_all
sample_data$d18Ow_mean_all <- pooled_data_stats$d18Ow_mean_all
sample_data$d18Ow_95CL_all <- pooled_data_stats$d18Ow_95CL_all

# ------------------------------------------------------------------------------------------
# Calculate mean temperature +/- 95%CL per horizon
# Group data by Horizon and calculate the mean and the pooled standard deviation for D47
pooled_horizon_stats <- data %>%
    subset(Type == "sample") %>%
    group_by(horizon) %>%
    summarise(
        N = n(),
        mean_D47 = mean(D47, na.rm = TRUE),
        D47_pooled_SD = sqrt(sum((D47 - mean(D47, na.rm = TRUE)) ** 2, na.rm = TRUE) / (n() - 1)),
        D47_se = D47_pooled_SD / sqrt(N),
        D47_CL95 = qt(0.975, N - 1) * D47_se,
        T47_horizon = sqrt(10 ** 6 / ((mean_D47 - 0.154) / 0.0391)) - 273.15, # Convert D47 to temperature using Anderson equation
        T47_max_horizon = sqrt(10 ** 6 / ((mean_D47 - D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 - 2 standard errors to temperature using Anderson equation
        T47_min_horizon = sqrt(10 ** 6 / ((mean_D47 + D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 + 2 standard errors to temperature using Anderson equation
        d18Ow_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_horizon) / -4.34, # Calculate d18Ow from mean d18Oc and mean T47_95CL_all
        d18Ow_max_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_max_horizon) / -4.34, # Calculate d18Ow from mean d18Oc and T47 - 2 standard errors
        d18Ow_min_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_min_horizon) / -4.34 # Calculate d18Ow from mean d18Oc and T47 + 2 standard errors
    ) %>%
    distinct(horizon, .keep_all = TRUE)

# Merge horizon statistics back into the original data
sample_data <- merge(
    sample_data,
    select(
        pooled_horizon_stats,
        horizon,
        T47_horizon,
        T47_max_horizon,
        T47_min_horizon,
        d18Ow_horizon,
        d18Ow_max_horizon,
        d18Ow_min_horizon
        ),
    by = "horizon"
)

# ------------------------------------------------------------------------------------------
# Plot temperatures per sample and add uncertainty bars per horizon
# Order the sample names first by horizon and then alphabetically
sample_data <- sample_data %>%
    arrange(horizon, Sample_name)

# Order the horizons in the specified order
sample_data$horizon <- factor(sample_data$horizon, levels = c("Zk13-6-1", "Zk13-6-4", "Zk13-21"))

# Plot temperatures per sample and add uncertainty bars per horizon
ggplot(sample_data,
    aes(
        x = factor(Sample_name, levels = unique(Sample_name)),
        y = T47_sample,
        color = factor(horizon, levels = c("Zk13-6-1", "Zk13-6-4", "Zk13-21")),
        fill = factor(horizon, levels = c("Zk13-6-1", "Zk13-6-4", "Zk13-21")),
        shape = Sample.description
    )
) +
geom_rect(aes(xmin = as.numeric(factor(Sample_name, levels = unique(Sample_name))) - 0.49, xmax = as.numeric(factor(Sample_name, levels = unique(Sample_name))) + 0.49, ymin = T47_min_horizon, ymax = T47_max_horizon), alpha = 0.5) +
geom_point() +
geom_errorbar(aes(ymin = T47_min_sample, ymax = T47_max_sample), width = 0.2, color = "black") +
labs(x = "Sample", y = "Temperature (°C)", title = "Temperature Estimates per Sample with 95% CL per Horizon", color = "horizon", fill = "horizon") +
scale_x_discrete(limits = unique(sample_data$Sample_name[order(sample_data$horizon)])) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1 / 1.5)

# ------------------------------------------------------------------------------------------
# Calculate mean temperature +/- 95%CL per season
# Group data by Season and calculate the mean and the pooled standard deviation for D47
pooled_season_stats <- data %>%
    subset(Type == "sample") %>%
    group_by(season) %>%
    summarise(
        N = n(),
        mean_D47 = mean(D47, na.rm = TRUE),
        D47_pooled_SD = sqrt(sum((D47 - mean(D47, na.rm = TRUE)) ** 2, na.rm = TRUE) / (n() - 1)),
        D47_se = D47_pooled_SD / sqrt(N),
        D47_CL95 = qt(0.975, N - 1) * D47_se,
        T47_season = sqrt(10 ** 6 / ((mean_D47 - 0.154) / 0.0391)) - 273.15, # Convert D47 to temperature using Anderson equation
        T47_max_season = sqrt(10 ** 6 / ((mean_D47 - D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 - 2 standard errors to temperature using Anderson equation
        T47_min_season = sqrt(10 ** 6 / ((mean_D47 + D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 + 2 standard errors to temperature using Anderson equation
        d18Ow_season = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_season) / -4.34, # Calculate d18Ow from mean d18Oc and mean T47_95CL_all
        d18Ow_max_season = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_max_season) / -4.34, # Calculate d18Ow from mean d18Oc and T47 - 2 standard errors
        d18Ow_min_season = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_min_season) / -4.34 # Calculate d18Ow from mean d18Oc and T47 + 2 standard errors
    )

# Merge season statistics back into the original data
sample_data <- merge(
    sample_data,
    select(
        pooled_season_stats,
        season,
        T47_season,
        T47_max_season,
        T47_min_season,
        d18Ow_season,
        d18Ow_max_season,
        d18Ow_min_season
        ),
    by = "season"
)

# ------------------------------------------------------------------------------------------
# Plot temperatures per sample and add uncertainty bars per season
# Order the sample names first by season and then alphabetically
sample_data <- sample_data %>%
    arrange(season, Sample_name)

# Order the seasons in the specified order
sample_data$season <- factor(sample_data$season, levels = c("first quarter", "third quarter", "no season"))

# Plot temperatures per sample and add uncertainty bars per season
ggplot(sample_data,
    aes(
        x = factor(Sample_name, levels = unique(Sample_name)),
        y = T47_sample,
        color = factor(season, levels = c("first quarter", "third quarter", "no season")),
        fill = factor(season, levels = c("first quarter", "third quarter", "no season")),
        shape = Sample.description
    )
) +
geom_rect(aes(xmin = as.numeric(factor(Sample_name, levels = unique(Sample_name))) - 0.49, xmax = as.numeric(factor(Sample_name, levels = unique(Sample_name))) + 0.49, ymin = T47_min_season, ymax = T47_max_season), alpha = 0.5) +
geom_point() +
geom_errorbar(aes(ymin = T47_min_sample, ymax = T47_max_sample), width = 0.2, color = "black") +
labs(x = "Sample", y = "Temperature (°C)", title = "Temperature Estimates per Sample with 95% CL per Season", color = "season", fill = "season") +
scale_x_discrete(limits = unique(sample_data$Sample_name[order(sample_data$season)])) +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1 / 1.5)

# ------------------------------------------------------------------------------------------
# Calculate mean temperature +/- 95%CL per season and horizon
# Group data by Season and Horizon and calculate the mean and the pooled standard deviation for D47
pooled_season_horizon_stats <- data %>%
    subset(Type == "sample") %>%
    group_by(season, horizon) %>%
    summarise(
        N = n(),
        mean_D47 = mean(D47, na.rm = TRUE),
        D47_pooled_SD = sqrt(sum((D47 - mean(D47, na.rm = TRUE)) ** 2, na.rm = TRUE) / (n() - 1)),
        D47_se = D47_pooled_SD / sqrt(N),
        D47_CL95 = qt(0.975, N - 1) * D47_se,
        T47_season_horizon = sqrt(10 ** 6 / ((mean_D47 - 0.154) / 0.0391)) - 273.15, # Convert D47 to temperature using Anderson equation
        T47_max_season_horizon = sqrt(10 ** 6 / ((mean_D47 - D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 - 2 standard errors to temperature using Anderson equation
        T47_min_season_horizon = sqrt(10 ** 6 / ((mean_D47 + D47_CL95 - 0.154) / 0.0391)) - 273.15, # Convert D47 + 2 standard errors to temperature using Anderson equation
        T47_95CL_season_horizon = (T47_max_season_horizon - T47_min_season_horizon) / 2, # Calculate 95% CL for temperature
        d18Ow_season_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_season_horizon) / -4.34, # Calculate d18Ow from mean d18Oc and mean T47_95CL_all
        d18Ow_max_season_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_max_season_horizon) / -4.34, # Calculate d18Ow from mean d18Oc and T47 - 2 standard errors
        d18Ow_min_season_horizon = (19.43 - 4.34 * mean(d18O, na.rm = TRUE) - T47_min_season_horizon) / -4.34, # Calculate d18Ow from mean d18Oc and T47 + 2 standard errors
        d18Ow_95CL_season_horizon = (d18Ow_max_season_horizon - d18Ow_min_season_horizon) / 2 # Calculate 95% CL for d18Ow
    )

# Merge season statistics back into the original data
sample_data <- merge(
    sample_data,
    select(
        pooled_season_horizon_stats,
        season,
        horizon,
        T47_season_horizon,
        T47_max_season_horizon,
        T47_min_season_horizon,
        T47_95CL_season_horizon,
        d18Ow_season_horizon,
        d18Ow_max_season_horizon,
        d18Ow_min_season_horizon,
        d18Ow_95CL_season_horizon
        ),
    by = c("season", "horizon")
)

# Export sample data
write.csv(sample_data, "sample_data.csv", row.names = FALSE)
