## Script to calculate individual PCB "personal" sampling rates
# for silicone wristbands.
# nd non-dominant hand
# d dominant hand

# Install packages
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("rollapply")

# Load libraries
{
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
  library(zoo)
}

# 3 volunteers, V1, V2, and V3
# Read data ---------------------------------------------------------------
{
  data.V1 <- read.csv("Data/Volunteer1.csv")
  data.V2 <- read.csv("Data/Volunteer2.csv")
  data.V3.1 <- read.csv("Data/Volunteer3.1.csv")
  data.V3.2 <- read.csv("Data/Volunteer3.2.csv")
  logKoa <- read.csv("Data/logKoa.csv")
  # ko from SamplingRates.R file
  ko <- read.csv("Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV1.csv")
  # Select only ko [m/d]
  ko <- ko[c(2,6)]
}

# Organize all dataset to have the same PCB congener list
# Get PCB names from each dataset
pcbs_V1 <- names(data.V1)[grep("^PCB", names(data.V1))]
pcbs_V2 <- names(data.V2)[grep("^PCB", names(data.V2))]
pcbs_V3.1 <- names(data.V3.1)[grep("^PCB", names(data.V3.1))]
pcbs_V3.2 <- names(data.V3.2)[grep("^PCB", names(data.V3.2))]
pcbs_ko <- ko$congener
pcbs_logKoa <- logKoa$congener

# Find common congeners across all
common_pcbs <- Reduce(intersect, list(pcbs_V1, pcbs_V2, pcbs_V3.1, pcbs_V3.2,
                                      pcbs_ko, pcbs_logKoa))
length(common_pcbs)

# Subset each dataset to include only the common PCBs
data.V1.pcbs <- data.V1[, c("time.day", "sample.code", common_pcbs)]
data.V2.pcbs <- data.V2[, c("time.day", "sample.code", common_pcbs)]
data.V3.1.pcbs <- data.V3.1[, c("time.day", "sample.code", common_pcbs)]
data.V3.2.pcbs <- data.V3.2[, c("time.day", "sample.code", common_pcbs)]

# Subset ko and logKoa data frames
ko.common <- ko[ko$congener %in% common_pcbs, ]
logKoa.common <- logKoa[logKoa$congener %in% common_pcbs, ]

# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa.common$congener,
  logKwb = 0.6156 * logKoa.common$logKoa + 2.161) # R2 = 0.96

# Calculate personal sampling rate V1 -------------------------------------
# Static WBs were used to calculate PCB concentration
# Triplicates for 4.3 days were deployed
# Effective volumes were calculated from the static ko, & Kws from above
# regression (logKws vs. logKoa)
# Both hands (d and nd)
{
  # Select WBs to calculate air concentration
  data.V1.1 <- data.V1.pcbs[1:3, ]
  # Average 3 WBs. NA values not included in the calculations
  data.V1.2 <- colMeans(data.V1.1[, 3:173], na.rm = TRUE)
  # Calculate air concentration in ng/m3
  # Use effective volume. Adult WBs
  Vwb <- data.V1$vol.WB[1]
  Awb <- data.V1$area.WB[1]
  # Calculate efective volume for static WBs
  veff_static.V1 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko.common$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.V1[1, 1]))
  # Compute concentration
  conc.V1 <- data.V1.2 / veff_static.V1
  # Calculate effective volume (Veff)
  subset_data <- data.V1.pcbs[4:13, 3:173]
  Veff.V1 <- t(apply(subset_data, 1, function(row) row / conc.V1))
  # Add metadata to Veff.V1 and change format
  Veff.V1 <- cbind(data.V1.pcbs[4:13, 2], data.V1.pcbs[4:13, 1], Veff.V1)
  # Transform to data.frame
  Veff.V1 <- as.data.frame(Veff.V1)
  # Add names to first 2 columns
  colnames(Veff.V1)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.V1[, 2:173] <- apply(Veff.V1[, 2:173], 2, as.numeric)
  # Select right, remove metadata
  Veff.V1.nd <- Veff.V1[1:5, 3:173]
  # Select time
  Veff.V1.nd.t <- Veff.V1[1:5, 2]
  # Select left, remove metadata
  Veff.V1.d <- Veff.V1[6:10, 3:173]
  # Select time
  Veff.V1.d.t <- Veff.V1[6:10, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.V1.nd <- matrix(nrow = length(Veff.V1.nd[1,]), ncol = 3)

for(i in 1:length(SR.V1.nd[, 1])) {
  if (length(unique(Veff.V1.nd[, i])) >= 3) {
    fit <- lm(Veff.V1.nd[, i] ~ 0 + Veff.V1.nd.t)
    SR.V1.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V1.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V1.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V1.nd[i, 1] <- 0
    SR.V1.nd[i, 2] <- 0
    SR.V1.nd[i, 3] <- 0
  }
}

SR.V1.nd <- data.frame(SR.V1.nd, group = "V1.nd")
colnames(SR.V1.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V1.nd)[0, ])
SR.V1.nd <- cbind(congener, SR.V1.nd)

# Convert R2 and p-value to numeric
SR.V1.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.V1.nd$`Sampling_Rate (m3/d)`)
SR.V1.nd$R2 <- as.numeric(SR.V1.nd$R2)
SR.V1.nd$p_value <- as.numeric(SR.V1.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V1.nd$R2 < 0.9 | SR.V1.nd$p_value > 0.05
SR.V1.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V1.nd$R2[mask] <- NA
SR.V1.nd$p_value[mask] <- NA
# Calculate ko from V1 nd
Awb.V1 <- data.V1$area.WB[4] # [m2] youth
SR.V1.nd$ko <- SR.V1.nd$`Sampling_Rate (m3/d)` / Awb.V1 # [m/d]

# Export results
write.csv(SR.V1.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V1.nd.csv", row.names = FALSE)

# Create matrix for sampling rate (SR)
SR.V1.d <- matrix(nrow = length(Veff.V1.d[1,]), ncol = 3)

for(i in 1:length(SR.V1.d[, 1])) {
  if (length(unique(Veff.V1.d[, i])) >= 3) {
    fit <- lm(Veff.V1.d[, i] ~ 0 + Veff.V1.d.t)
    SR.V1.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V1.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V1.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V1.d[i, 1] <- 0
    SR.V1.d[i, 2] <- 0
    SR.V1.d[i, 3] <- 0
  }
}

SR.V1.d <- data.frame(SR.V1.d, group = "V1.d")
colnames(SR.V1.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V1.d)[0, ])
SR.V1.d <- cbind(congener, SR.V1.d)

# Convert R2 and p-value to numeric
SR.V1.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.V1.d$`Sampling_Rate (m3/d)`)
SR.V1.d$R2 <- as.numeric(SR.V1.d$R2)
SR.V1.d$p_value <- as.numeric(SR.V1.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V1.d$R2 < 0.9 | SR.V1.d$p_value > 0.05
SR.V1.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V1.d$R2[mask] <- NA
SR.V1.d$p_value[mask] <- NA
# Calculate ko from V1 d
# Calculate ko from V1 nd
Awb.V1 <- data.V1$area.WB[4] # [m2] youth
SR.V1.d$ko <- SR.V1.d$`Sampling_Rate (m3/d)` / Awb.V1 # [m/d]

# Export results
write.csv(SR.V1.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V1.d.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.V1.nd$congener <- factor(SR.V1.nd$congener,
                            levels = unique(SR.V1.nd$congener))
# Plot with legend
ggplot(SR.V1.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Organize PCB names
SR.V1.d$congener <- factor(SR.V1.d$congener,
                               levels = unique(SR.V1.d$congener))
# Plot with legend
ggplot(SR.V1.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# V1 SR vs logKoa regression ------------------------------------------
# (1) Average both d and nd
sr.ave.V1 <- as.data.frame(rowMeans(cbind(SR.V1.d$`Sampling_Rate (m3/d)`, 
                                SR.V1.nd$`Sampling_Rate (m3/d)`), 
                          na.rm = TRUE))

sr.ave.V1$logkoa <- logKoa.common$logKoa
colnames(sr.ave.V1) <- c('ave_sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V1.1 <- lm(log(sr.ave.V1$ave_sr) ~ sr.ave.V1$logKoa)

# Get the coefficients
a <- exp(coef(model.V1.1)[1])  # exponentiate the intercept
b <- coef(model.V1.1)[2]       # coefficient for logKoa
r2 <- summary(model.V1.1)$r.squared

# plot
p.sr.V1.koa.1 <- ggplot(sr.ave.V1, aes(x = logKoa, y = ave_sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7, y = 15.7,
           label = paste("Ave. Vol. 1 (d & nd)"),size = 4) +
  annotate("text", x = 7.5, y = 15,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 4) + 
  annotate("text", x = 6.6, y = 14.2,
           label = paste("R² = ", round(r2, 2)), size = 4) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.V1.koa.1

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V1_logKoa1.png", plot = p.sr.V1.koa.1,
       width = 6, height = 6, dpi = 500)

# (2) Individual values
# Create a long dataframe combining SR.V1.d and SR.V1.nd
sr.long.V1 <- data.frame(
  sr = c(SR.V1.d$`Sampling_Rate (m3/d)`, SR.V1.nd$`Sampling_Rate (m3/d)`),
  logKoa = rep(logKoa.common$logKoa, 2)  # Repeat logKoa values for both d and nd
)

# Remove any NA values
sr.long.V1 <- na.omit(sr.long.V1)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V1.2 <- lm(log(sr.long.V1$sr) ~ sr.long.V1$logKoa)

# Get the coefficients
a <- exp(coef(model.V1.2)[1])  # exponentiate the intercept
b <- coef(model.V1.2)[2]       # coefficient for logKoa
r2 <- summary(model.V1.2)$r.squared

# Plot
p.sr.V1.koa.2 <- ggplot(sr.long.V1, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.long.V1$logKoa) + 0.6, y = max(sr.long.V1$sr) * 1.2,
           label = paste("Vol. 1 (d & nd)"),size = 5) +
  annotate("text", x = min(sr.long.V1$logKoa) + 1.5, y = max(sr.long.V1$sr) * 1.15,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.long.V1$logKoa) + 0.35, y = max(sr.long.V1$sr) * 1.1,
           label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V1.koa.2

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V1_logKoa2.png", plot = p.sr.V1.koa.2,
       width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate V2 -------------------------------------
# WBs were used to calculate PCB concentration
# Dominant hand (d)
# triplicates for 4.27 days were deployed
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.V2.1 <- data.V2.pcbs[1:3,]
  # Average 3 WBs
  data.V2.2 <- colMeans(data.V2.1[, 3:173])
  # Calculate air concentration in ng/m3
  # Use effective volume. Adult WBs
  Vwb <- data.V1$vol.WB[1]
  Awb <- data.V1$area.WB[1]
  veff_stat.V2 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko.common$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.V2[1, 1]))
  # Compute concentration
  conc.V2 <- data.V2.2 / veff_stat.V2
  # Calculate effective volume (Veff)
  subset_data <- data.V2.pcbs[4:8, 3:173]
  Veff.V2 <- t(apply(subset_data, 1, function(row) row / conc.V2))
  # Add metadata to Veff.V2 and change format
  Veff.V2 <- cbind(data.V2.pcbs[4:8, 2], data.V2.pcbs[4:8, 1], Veff.V2)
  # Transform to data.frame
  Veff.V2 <- as.data.frame(Veff.V2)
  # Add names to first 2 columns
  colnames(Veff.V2)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.V2[, 2:173] <- apply(Veff.V2[, 2:173], 2, as.numeric)
  # Select right, remove metadata
  Veff.V2.d <- Veff.V2[1:5, 3:173]
  # Select time
  Veff.V2.d.t <- Veff.V2[1:5, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.V2.d <- matrix(nrow = length(Veff.V2.d[1,]), ncol = 3)

for(i in 1:length(SR.V2.d[, 1])) {
  if (length(unique(Veff.V2.d[, i])) >= 3) {
    fit <- lm(Veff.V2.d[, i] ~ 0 + Veff.V2.d.t)
    SR.V2.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V2.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V2.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V2.d[i, 1] <- 0
    SR.V2.d[i, 2] <- 0
    SR.V2.d[i, 3] <- 0
  }
}

SR.V2.d <- data.frame(SR.V2.d, group = "V2.d")
colnames(SR.V2.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V2.d)[0, ])
SR.V2.d <- cbind(congener, SR.V2.d)

# Convert R2 and p-value to numeric
SR.V2.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.V2.d$`Sampling_Rate (m3/d)`)
SR.V2.d$R2 <- as.numeric(SR.V2.d$R2)
SR.V2.d$p_value <- as.numeric(SR.V2.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V2.d$R2 < 0.9 | SR.V2.d$p_value > 0.05
SR.V2.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V2.d$R2[mask] <- NA
SR.V2.d$p_value[mask] <- NA
# Calculate ko from V2 d
Awb.V2 <- data.V2$area.WB[1] # [m2] adult
SR.V2.d$ko <- SR.V2.d$`Sampling_Rate (m3/d)` / Awb.V2 # [m/d]

# Export results
write.csv(SR.V2.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V2.d.csv", row.names = FALSE)

# Plot
# Organize PCB names
SR.V2.d$congener <- factor(SR.V2.d$congener,
                               levels = unique(SR.V2.d$congener))

ggplot(SR.V2.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# V2 SR vs logKoa regression ---------------------------------------------
# Create a long dataframe combining SR.V2.d
sr.V2 <- data.frame(
  sr = SR.V2.d$`Sampling_Rate (m3/d)`,
  logKoa = logKoa.common$logKoa)

# Remove any NA values
sr.V2 <- na.omit(sr.V2)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V2 <- lm(log(sr.V2$sr) ~ sr.V2$logKoa)

# Get the coefficients
a <- exp(coef(model.V2)[1])  # exponentiate the intercept
b <- coef(model.V2)[2]       # coefficient for logKoa
r2 <- summary(model.V2)$r.squared

# Plot
p.sr.V2.koa <- ggplot(sr.V2, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 6.6, y = 5, label = paste("Vol. 2 (d)"),size = 5) +
  annotate("text", x = 7.67, y = 4.7,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = 6.58, y = 4.4, label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V2.koa

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V2_logKoa.png", plot = p.sr.V2.koa,
       width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate V3.1 -----------------------------------
# WBs were used to calculate PCB concentration
# Non-dominant hand (nd)
# Concentrations were calculated for each sampling day 
{
  # Constants
  Vwb <- data.V1$vol.WB[1]
  Awb <- data.V1$area.WB[1]
  # Prepare data
  data.V3.1.1 <- data.V3.1.pcbs[1:6, 3:173]
  Kwb_val <- 10^logKwb$logKwb
  ko_val <- ko.common$ko
  # Initialize lists to store results
  veff_stat_list <- list()
  conc_list <- list()
  Veff_list <- list()
  for (i in 1:6) {
    # Effective volume calculation
    veff_stat <- Kwb_val * Vwb * (1 - exp(-ko_val * Awb / Vwb / Kwb_val * data.V3.1[i, 1]))
    veff_stat_list[[i]] <- veff_stat
    # Air concentration calculation
    conc <- data.V3.1.1[i, ] / veff_stat
    conc_list[[i]] <- conc
    # Effective volume (Veff)
    Veff <- data.V3.1.pcbs[i + 6, 3:173] / conc
    Veff_list[[i]] <- Veff
  }
  
  # Combine Veff results
  Veff.V3.1 <- do.call(rbind, Veff_list)
  Veff.V3.1 <- cbind(data.V3.1[7:12, 1:2], Veff.V3.1)
  Veff.V3.1 <- as.data.frame(Veff.V3.1)
  # Select weeks
  Veff.V3.1st.nd <- Veff.V3.1[1:3, 3:173]
  Veff.V3.1st.nd.t <- Veff.V3.1[1:3, 1]
  Veff.V3.2nd.nd <- Veff.V3.1[4:6, 3:173]
  Veff.V3.2nd.nd.t <- Veff.V3.1[4:6, 1]
}

# Calculate sampling rate (SR) for 1st and 2nd weeks (m3/d)
# Create matrix for sampling rate (SR)
SR.V3.1st.nd <- matrix(nrow = length(Veff.V3.1st.nd[1,]), ncol = 3)

for(i in 1:length(SR.V3.1st.nd[, 1])) {
  if (sum(!is.na(Veff.V3.1st.nd[, i ]) & !is.infinite(Veff.V3.1st.nd[, i])) == 3) {
    fit <- lm(Veff.V3.1st.nd[, i] ~ 0 + Veff.V3.1st.nd.t)
    SR.V3.1st.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V3.1st.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V3.1st.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V3.1st.nd[i, 1] <- 0
    SR.V3.1st.nd[i, 2] <- 0
    SR.V3.1st.nd[i, 3] <- 0
  }
}

SR.V3.1st.nd <- data.frame(SR.V3.1st.nd, group = "V3.1st.nd")
colnames(SR.V3.1st.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V3.1st.nd)[0, ])
SR.V3.1st.nd <- cbind(congener, SR.V3.1st.nd)

# Convert R2 and p-value to numeric
SR.V3.1st.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.V3.1st.nd$`Sampling_Rate (m3/d)`)
SR.V3.1st.nd$R2 <- as.numeric(SR.V3.1st.nd$R2)
SR.V3.1st.nd$p_value <- as.numeric(SR.V3.1st.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V3.1st.nd$R2 < 0.9 | SR.V3.1st.nd$p_value > 0.05
SR.V3.1st.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V3.1st.nd$R2[mask] <- NA
SR.V3.1st.nd$p_value[mask] <- NA
# Calculate ko from V3.1
Awb.V3.1 <- data.V3.1$area.WB[1] # [m2] adult
SR.V3.1st.nd$ko <- SR.V3.1st.nd$`Sampling_Rate (m3/d)` / Awb.V3.1 # [m/d]

# Export results
write.csv(SR.V3.1st.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V3.1st.nd.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.V3.1st.nd$congener <- factor(SR.V3.1st.nd$congener,
                            levels = unique(SR.V3.1st.nd$congener))

ggplot(SR.V3.1st.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Create matrix for sampling rate (SR)
SR.V3.2nd.nd <- matrix(nrow = length(Veff.V3.2nd.nd[1,]), ncol = 3)

for(i in 1:length(SR.V3.2nd.nd[, 1])) {
  if (sum(!is.na(Veff.V3.2nd.nd[, i ]) & !is.infinite(Veff.V3.2nd.nd[, i])) == 3) {
    fit <- lm(Veff.V3.2nd.nd[, i] ~ 0 + Veff.V3.2nd.nd.t)
    SR.V3.2nd.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V3.2nd.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V3.2nd.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V3.2nd.nd[i, 1] <- 0
    SR.V3.2nd.nd[i, 2] <- 0
    SR.V3.2nd.nd[i, 3] <- 0
  }
}

SR.V3.2nd.nd <- data.frame(SR.V3.2nd.nd, group = "V3.2nd.nd")
colnames(SR.V3.2nd.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V3.2nd.nd)[0, ])
SR.V3.2nd.nd <- cbind(congener, SR.V3.2nd.nd)

# Convert R2 and p-value to numeric
SR.V3.2nd.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.V3.2nd.nd$`Sampling_Rate (m3/d)`)
SR.V3.2nd.nd$R2 <- as.numeric(SR.V3.2nd.nd$R2)
SR.V3.2nd.nd$p_value <- as.numeric(SR.V3.2nd.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V3.2nd.nd$R2 < 0.9 | SR.V3.2nd.nd$p_value > 0.05
SR.V3.2nd.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V3.2nd.nd$R2[mask] <- NA
SR.V3.2nd.nd$p_value[mask] <- NA
# Calculate ko from V3.1
Awb.V3.2nd <- data.V3.1$area.WB[1] # [m2] adult
SR.V3.2nd.nd$ko <- SR.V3.2nd.nd$`Sampling_Rate (m3/d)` / Awb.V3.2nd # [m/d]

# Export results
write.csv(SR.V3.2nd.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V3.2nd.nd.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.V3.2nd.nd$congener <- factor(SR.V3.2nd.nd$congener,
                              levels = unique(SR.V3.2nd.nd$congener))

ggplot(SR.V3.2nd.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# V3 SR vs logKoa regression 1 ----------------------------------------------
# (1) Average both nds
sr.ave.V3 <- as.data.frame(rowMeans(cbind(SR.V3.1st.nd$`Sampling_Rate (m3/d)`, 
                                              SR.V3.2nd.nd$`Sampling_Rate (m3/d)`), 
                                        na.rm = TRUE))

sr.ave.V3$logkoa <- logKoa.common$logKoa
colnames(sr.ave.V3) <- c('ave_sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.1 <- lm(log(sr.ave.V3$ave_sr) ~ sr.ave.V3$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.1)[1])  # exponentiate the intercept
b <- coef(model.V3.1)[2]       # coefficient for logKoa
r2 <- summary(model.V3.1)$r.squared

# plot
p.sr.V3.koa.1 <- ggplot(sr.ave.V3, aes(x = logKoa, y = ave_sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.6, y = 7,
           label = paste("Ave. Vol. 3 (nd 1st & 2nd weeks)"),size = 5) +
  annotate("text", x = 7.55, y = 6.7,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 5) + 
  annotate("text", x = 6.75, y = 6.4,
           label = paste("R² = ", round(r2, 2)), size = 5) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Average Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V3.koa.1

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.nd.1.png",
       plot = p.sr.V3.koa.1, width = 6, height = 6, dpi = 500)

# (2) Individual values
# Create a long dataframe combining SR.V3.1st. nd and SR.V3.2nd.nd
sr.long.V3 <- data.frame(
  sr = c(SR.V3.1st.nd$`Sampling_Rate (m3/d)`, SR.V3.2nd.nd$`Sampling_Rate (m3/d)`),
  logKoa = rep(logKoa.common$logKoa, 2)  # Repeat logKoa values for both d and nd
)

# Remove any NA values
sr.long.V3 <- na.omit(sr.long.V3)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.2 <- lm(log(sr.long.V3$sr) ~ sr.long.V3$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.2)[1])  # exponentiate the intercept
b <- coef(model.V3.2)[2]       # coefficient for logKoa
r2 <- summary(model.V3.2)$r.squared

# Plot
p.sr.V3.koa.2 <- ggplot(sr.long.V3, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.8, y = 7,
           label = paste("Vol. 3 (nd 1st & 2nd weeks)"),size = 5) +
  annotate("text", x = min(sr.long.V3$logKoa) + 1.45, y = max(sr.long.V3$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.long.V3$logKoa) + 0.35, y = max(sr.long.V3$sr) * 1.13,
           label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.V3.koa.2

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.nd.2.png",
       plot = p.sr.V3.koa.2, width = 6, height = 6, dpi = 500)

# (3) 1st week only
sr.long.V3.1st <- data.frame(
  sr = c(SR.V3.1st.nd$`Sampling_Rate (m3/d)`),
  logKoa = logKoa.common$logKoa)

# Remove any NA values
sr.long.V3.1st <- na.omit(sr.long.V3.1st)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.3 <- lm(log(sr.long.V3.1st$sr) ~ sr.long.V3.1st$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.3)[1])  # exponentiate the intercept
b <- coef(model.V3.3)[2]       # coefficient for logKoa
r2 <- summary(model.V3.3)$r.squared

# Plot
p.sr.V3.koa.3 <- ggplot(sr.long.V3.1st, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  xlim(7, 11) +
  annotate("text", x = 7.6, y = 7,
           label = paste("Vol. 3 (nd 1st week)"),size = 5) +
  annotate("text", x = 8, y = 6.7,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 5) + 
  annotate("text", x = 7.25, y = 6.4,
           label = paste("R² = ", round(r2, 2)), size = 5) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V3.koa.3

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.nd.1st.week.png",
       plot = p.sr.V3.koa.3, width = 6, height = 6, dpi = 500)

# (4) 2nd week only
sr.long.V3.2nd <- data.frame(
  sr = c(SR.V3.2nd.nd$`Sampling_Rate (m3/d)`),
  logKoa = logKoa.common$logKoa)

# Remove any NA values
sr.long.V3.2nd <- na.omit(sr.long.V3.2nd)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.4 <- lm(log(sr.long.V3.2nd$sr) ~ sr.long.V3.2nd$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.4)[1])  # exponentiate the intercept
b <- coef(model.V3.4)[2]       # coefficient for logKoa
r2 <- summary(model.V3.4)$r.squared

# Plot
p.sr.V3.koa.4 <- ggplot(sr.long.V3.2nd, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.4, y = 4.5,
           label = paste("Vol. 3 (nd 2nd week)"),size = 5) +
  annotate("text", x = 7.9, y = 4.3,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 5) + 
  annotate("text", x = 6.95, y = 4.1,
           label = paste("R² = ", round(r2, 2)), size = 5) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V3.koa.4

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.nd.2nd.week.png",
       plot = p.sr.V3.koa.4, width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate V3.2 -----------------------------------
# Need to check the code!!
# 3 WBs were not wiped and 3 WBs were wiped
# Dominant hand used here
# WBs were used to calculate PCB concentration
# Concentrations were calculated for each sampling day 
# sampling rate of 0.5 m3/d was used for static WBs

# Constants
Vwb <- data.V3.2$vol.WB[1]
Awb <- data.V3.2$area.WB[1]
Kwb_val <- 10^logKwb$logKwb  # Assuming logKwb is available for all congeners
ko_val <- ko.common$ko  # Assuming ko values are available for all congeners

# Initialize lists to store results for each PCB
SR.V3.2_results <- list()

# Loop over each PCB (assuming PCB names are in `logKoa.common$congener`)
for (congener in logKoa.common$congener) {
  
  # Extract values for current PCB
  PCB_vals <- data.V3.2.pcbs[[congener]]
  logKoa_val <- logKoa.common$logKoa[logKoa.common$congener == congener]
  ko_val_for_PCB <- ko_val[ko.common$congener == congener]
  
  # Time values for the first three rows
  time_vals <- data.V3.2$time.day[1:3]
  
  # Calculate Kwb
  Kwb_val_for_PCB <- 10^(0.6156 * logKoa_val + 2.161)
  
  # Calculate veff_stat for 3 static rows
  veff_stat <- sapply(time_vals, function(time) {
    Kwb_val_for_PCB * Vwb * (1 - exp(-ko_val_for_PCB * Awb * time / Vwb / Kwb_val_for_PCB))
  })
  
  # Calculate concentration using static rows (1 to 3)
  conc <- PCB_vals[1:3] / veff_stat
  
  # Calculate Veff for rows 4 to 9
  Veff <- sapply(4:9, function(row) {
    cycle_index <- ((row - 4) %% 3) + 1
    PCB_vals[row] / conc[cycle_index]
  })
  
  # Create data frame for Veff results
  Veff_results <- data.frame(
    Row = 4:9,
    Time = data.V3.2$time.day[4:9],
    Veff = Veff
  )
  
  # Split data into `nw` and `w` for regression
  Veff_nw <- Veff_results[1:3, ]
  Veff_w <- Veff_results[4:6, ]
  
  # Regression function
  get_regression <- function(values, times) {
    if (sum(!is.na(values)) == 3) {  # Ensure no NA or infinite values
      fit <- lm(values ~ 0 + times)
      SR <- format(signif(coef(summary(fit))[1, "Estimate"], digits = 3))
      R2 <- format(signif(summary(fit)$adj.r.squared, digits = 3))
      pval <- format(signif(coef(summary(fit))[1, "Pr(>|t|)"], digits = 3))
    } else {
      SR <- R2 <- pval <- "0"
    }
    return(c(SR, R2, pval))
  }
  
  # Run regressions for nw and w
  SR_nw <- get_regression(Veff_nw$Veff, Veff_nw$Time)
  SR_w <- get_regression(Veff_w$Veff, Veff_w$Time)
  
  # Store regression results in a data frame
  SR.V3.2_results[[congener]] <- data.frame(
    Sampling_Rate = c(SR_nw[1], SR_w[1]),
    R2 = c(SR_nw[2], SR_w[2]),
    p_value = c(SR_nw[3], SR_w[3]),
    group = c(paste(congener, "nw", sep = "_"), paste(congener, "w", sep = "_"))
  )
}

# Combine all results for all congeners
SR_V3.2 <- do.call(rbind, SR.V3.2_results)
# Change column name
colnames(SR_V3.2)[colnames(SR_V3.2) == "Sampling_Rate"] <- "Sampling_Rate (m3/d)"

# Convert columns to numeric
SR_V3.2$'Sampling_Rate (m3/d)' <- as.numeric(SR_V3.2$`Sampling_Rate (m3/d)`)
SR_V3.2$R2 <- as.numeric(SR_V3.2$R2)
SR_V3.2$p_value <- as.numeric(SR_V3.2$p_value)

# Define masks
mask_filter <- SR_V3.2$R2 < 0.9 | SR_V3.2$p_value > 0.05

# Apply NA to filtered rows
SR_V3.2[mask_filter, c("Sampling_Rate (m3/d)", "R2", "p_value")] <- NA

# Define areas
Awb.V3.nw <- data.V3.2$area.WB[3]  # NW
Awb.V3.w <- data.V3.2$area.WB[4]   # W

# Assign area based on group (using grepl for pattern matching)
SR_V3.2$area <- ifelse(grepl("_nw$", SR_V3.2$group), Awb.V3.nw,
                       ifelse(grepl("_w$", SR_V3.2$group), Awb.V3.w, NA))

# Calculate ko
SR_V3.2$ko <- SR_V3.2$Sampling_Rate / SR_V3.2$area

# Optional: drop the `area` column if not needed
SR_V3.2$area <- NULL

# Subset for nw and w
SR_V3.2_nw <- subset(SR_V3.2, grepl("_nw$", group))
SR_V3.2_w  <- subset(SR_V3.2, grepl("_w$", group))

SR_V3.2_nw$congener <- common_pcbs
SR_V3.2_w$congener <- common_pcbs
rownames(SR_V3.2_nw) <- NULL
rownames(SR_V3.2_w) <- NULL
SR_V3.2_nw <- SR_V3.2_nw[, c("congener", setdiff(names(SR_V3.2_nw), "congener"))]
SR_V3.2_w <- SR_V3.2_w[, c("congener", setdiff(names(SR_V3.2_w), "congener"))]
SR_V3.2_nw$group <- "V3.nw.d"
SR_V3.2_w$group <- "V3.w.d"

# Save to separate CSV files
write.csv(SR_V3.2_nw, "Output/Data/csv/SamplingRates/Personal/SR.V3.2_nw.csv", row.names = FALSE)
write.csv(SR_V3.2_w,  "Output/Data/csv/SamplingRates/Personal/SR.V3.2_w.csv",  row.names = FALSE)

# Plot
# (1) nw
# Organize PCB names
SR_V3.2_nw$congener <- factor(SR_V3.2_nw$congener,
                              levels = unique(SR_V3.2_nw$congener))

ggplot(SR_V3.2_nw, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# (2) w
SR_V3.2_w$congener <- factor(SR_V3.2_w$congener,
                              levels = unique(SR_V3.2_w$congener))

ggplot(SR_V3.2_w, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# V3 SR vs logKoa regression 2 ------------------------------------------
# Check difference btw w and nw, both d
# (1) nw
sr.V3.nw <- as.data.frame(SR_V3.2_nw$`Sampling_Rate`, na.rm = TRUE)

sr.V3.nw$logkoa <- logKoa.common$logKoa
colnames(sr.V3.nw) <- c('sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.nw <- lm(log(sr.V3.nw$sr) ~ sr.V3.nw$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.nw)[1])  # exponentiate the intercept
b <- coef(model.V3.nw)[2]       # coefficient for logKoa
r2 <- summary(model.V3.nw)$r.squared

# plot
p.sr.V3.koa.nw <- ggplot(sr.V3.nw, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 6.9, y = 5.4, label = paste("Vol. 3 (non-wiped)"),
           size = 5) +
  annotate("text", x = 7.6, y = 5.1,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 5) + 
  annotate("text", x = 6.52, y = 4.8,
           label = paste("R² = ", round(r2, 2)), size = 5) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.V3.koa.nw

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.nw.png",
       plot = p.sr.V3.koa.nw, width = 6, height = 6, dpi = 500)

# (2) w
sr.V3.w <- as.data.frame(SR_V3.2_w$`Sampling_Rate`, na.rm = TRUE)

sr.V3.w$logkoa <- logKoa.common$logKoa
colnames(sr.V3.w) <- c('sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.V3.w <- lm(log(sr.V3.w$sr) ~ sr.V3.w$logKoa)

# Get the coefficients
a <- exp(coef(model.V3.w)[1])  # exponentiate the intercept
b <- coef(model.V3.w)[2]       # coefficient for logKoa
r2 <- summary(model.V3.w)$r.squared

# plot
p.sr.V3.koa.w <- ggplot(sr.V3.w, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 6.7, y = 6.2, label = paste("Vol. 3 (wiped)")) +
  annotate("text", x = 7.5, y = 6,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 4) + 
  annotate("text", x = 6.6, y = 5.8,
           label = paste("R² = ", round(r2, 2)), size = 4) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.V3.koa.w

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/V3_logKoa.w.png",
       plot = p.sr.V3.koa.w, width = 6, height = 6, dpi = 500)

# Plot individual congeners -----------------------------------------------
# Combine plot
# Only volunteers 1 and 2 with five samples
# Combine data, padding shorter dataset with NA
max_len <- max(length(Veff.V1.d.t), length(Veff.V1.nd.t),
               length(Veff.V2.d.t))

# PCBs 18+30
# Plot included in the manuscript
up.V1.d <- data.frame(
  time = c(Veff.V1.d.t, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  veff = c(Veff.V1.d$PCB18.30, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.V1.nd <- data.frame(
  time = c(Veff.V1.nd.t, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  veff = c(Veff.V1.nd$PCB18.30, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.V2.d <- data.frame(
  time = c(Veff.V2.d.t, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  veff = c(Veff.V2.d$PCB18.30, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  group = rep("Vol. 2 d")
)

combined_data <- rbind(up.V1.d, up.V1.nd, up.V2.d)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(veff ~ 0 + time))[1]) # [m3/d]

# Define colors for groups
group_colors <- c("Vol. 1 d" = "#E65C00",
                  "Vol. 1 nd" = "#0072B2",
                  "Vol. 2 d" = "#009E73")

# Plot with defined colors
plot.18.30 <- ggplot(combined_data, aes(x = time * 24, y = veff,
                                        color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  theme_bw() +
  xlim(0, 45) +
  ylim(0, 2.0) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCBs 18+30 (m"^"3"*")"))) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold", size = 28),
        axis.title.y = element_text(face = "bold", size = 28),
        axis.text.x = element_text(face = "bold", size = 28),
        axis.title.x = element_text(face = "bold", size = 28),
        legend.position = "none",  # Remove the existing legend
        aspect.ratio = 1.5)

# Adding color symbols before text with increased spacing
plot.18.30 <- plot.18.30 + 
  geom_point(aes(x = 1, y = 2.0), size = 6, shape = 21,
             fill = group_colors["Vol. 1 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 2.0,
           label = bquote(Vol. ~ "1 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 1.85), size = 6, shape = 21,
             fill = group_colors["Vol. 1 nd"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 1.85,
           label = bquote(Vol. ~ "1 nd" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 nd"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 1.7), size = 6, shape = 21,
             fill = group_colors["Vol. 2 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 1.7,
           label = bquote(Vol. ~ "2 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 2 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 35, y = 0.1,
           label = "Study 2", hjust = 0, size = 10, color = "black") +
  annotate("text", x = Inf, y = Inf,
           label = "(b)", hjust = 1.1, vjust = 1.5, 
           size = 10, color = "black")

# See plot
plot.18.30

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB18.30VoluntSamplingRatesv5.png",
       plot = plot.18.30, width = 8, height = 10, dpi = 1300)

# PCB 52
# Plot included in the manuscript
up.V1.d <- data.frame(
  time = c(Veff.V1.d.t, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  veff = c(Veff.V1.d$PCB52, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.V1.nd <- data.frame(
  time = c(Veff.V1.nd.t, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  veff = c(Veff.V1.nd$PCB52, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.V2.d <- data.frame(
  time = c(Veff.V2.d.t, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  veff = c(Veff.V2.d$PCB52, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  group = rep("Vol. 2 d")
)

combined_data <- rbind(up.V1.d, up.V1.nd, up.V2.d)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(veff ~ 0 + time))[1]) # [m3/d]

# Define colors for groups
group_colors <- c("Vol. 1 d" = "#E65C00",
                  "Vol. 1 nd" = "#0072B2",
                  "Vol. 2 d" = "#009E73")

# Plot with defined colors
plot.52 <- ggplot(combined_data, aes(x = time * 24, y = veff,
                                     color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  theme_bw() +
  xlim(0, 45) +
  ylim(0, 3.1) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 52 (m"^"3"*")"))) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold", size = 28),
        axis.title.y = element_text(face = "bold", size = 28),
        axis.text.x = element_text(face = "bold", size = 28),
        axis.title.x = element_text(face = "bold", size = 28),
        legend.position = "none",
        aspect.ratio = 1.5)

# Adding color symbols before text with increased spacing
plot.52 <- plot.52 + 
  geom_point(aes(x = 1, y = 3.1), size = 6, shape = 21,
             fill = group_colors["Vol. 1 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 3.1,
           label = bquote(Vol. ~ "1 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 2.9), size = 6, shape = 21,
             fill = group_colors["Vol. 1 nd"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 2.9,
           label = bquote(Vol. ~ "1 nd" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 nd"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 2.7), size = 6, shape = 21,
             fill = group_colors["Vol. 2 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 2.7,
           label = bquote(Vol. ~ "2 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 2 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 35, y = 0.1,
           label = "Study 2", hjust = 0, size = 10, color = "black") +
  annotate("text", x = Inf, y = Inf, label = "(c)", hjust = 1.1, vjust = 1.5, 
           size = 10, color = "black")

# See plot
plot.52

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB52VoluntSamplingRatesV5.png",
       plot = plot.52, width = 8, height = 10, dpi = 1300)

# PCB 118
up.V1.d <- data.frame(
  time = c(Veff.V1.d.t, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  veff = c(Veff.V1.d$PCB118, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.V1.nd <- data.frame(
  time = c(Veff.V1.nd.t, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  veff = c(Veff.V1.nd$PCB118, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.V2.d <- data.frame(
  time = c(Veff.V2.d.t, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  veff = c(Veff.V2.d$PCB118, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  group = rep("Vol. 2 d")
)

combined_data <- rbind(up.V1.d, up.V1.nd, up.V2.d)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(veff ~ 0 + time))[1]) # [m3/d]

# Define colors for groups
group_colors <- c("Vol. 1 d" = "#E65C00",
                  "Vol. 1 nd" = "#0072B2",
                  "Vol. 2 d" = "#009E73")

# Plot with defined colors
plot.118 <- ggplot(combined_data, aes(x = time * 24, y = veff,
                                     color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  theme_bw() +
  xlim(0, 45) +
  ylim(0, 10) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 118 (m"^"3"*")"))) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
        legend.position = "none",
        aspect.ratio = 1.5)

# Adding color symbols before text with increased spacing
plot.118 <- plot.118 + 
  geom_point(aes(x = 1, y = 10.0), size = 6, shape = 21,
             fill = group_colors["Vol. 1 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 10.0,
           label = bquote(Vol. ~ "1 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 9.25), size = 6, shape = 21,
             fill = group_colors["Vol. 1 nd"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 9.25,
           label = bquote(Vol. ~ "1 nd" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 nd"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 8.5), size = 6, shape = 21,
             fill = group_colors["Vol. 2 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 8.5,
           label = bquote(Vol. ~ "2 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 2 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 35, y = 0.1,
           label = "Study 2", hjust = 0, size = 10, color = "black")

# See plot
plot.118

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB118VoluntSamplingRatesV2.png",
       plot = plot.118, width = 8, height = 10, dpi = 1300)

# PCB 187
# Plot included in the manuscript
up.V1.d <- data.frame(
  time = c(Veff.V1.d.t, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  veff = c(Veff.V1.d$PCB187, rep(NA, max_len - length(Veff.V1.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.V1.nd <- data.frame(
  time = c(Veff.V1.nd.t, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  veff = c(Veff.V1.nd$PCB187, rep(NA, max_len - length(Veff.V1.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.V2.d <- data.frame(
  time = c(Veff.V2.d.t, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  veff = c(Veff.V2.d$PCB187, rep(NA, max_len - length(Veff.V2.d[, 1]))),
  group = rep("Vol. 2 d")
)

combined_data <- rbind(up.V1.d, up.V1.nd, up.V2.d)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(veff ~ 0 + time))[1]) # [m3/d]

# Define colors for groups
group_colors <- c("Vol. 1 d" = "#E65C00",
                  "Vol. 1 nd" = "#0072B2",
                  "Vol. 2 d" = "#009E73")

# Plot with defined colors
plot.187 <- ggplot(combined_data, aes(x = time * 24, y = veff,
                                      color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  theme_bw() +
  xlim(0, 45) +
  ylim(0, 16) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 187 (m"^"3"*")"))) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold", size = 28),
        axis.title.y = element_text(face = "bold", size = 28),
        axis.text.x = element_text(face = "bold", size = 28),
        axis.title.x = element_text(face = "bold", size = 28),
        legend.position = "none",
        aspect.ratio = 1.5)

# Adding color symbols before text with increased spacing
plot.187 <- plot.187 + 
  geom_point(aes(x = 1, y = 16.0), size = 6, shape = 21,
             fill = group_colors["Vol. 1 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 16.0,
           label = bquote(Vol. ~ "1 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 15.0), size = 6, shape = 21,
             fill = group_colors["Vol. 1 nd"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 15.0,
           label = bquote(Vol. ~ "1 nd" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 1 nd"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  geom_point(aes(x = 1, y = 14.0), size = 6, shape = 21,
             fill = group_colors["Vol. 2 d"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 2, y = 14.0,
           label = bquote(Vol. ~ "2 d" ~ "=" ~ .(round(slopes$slope[slopes$group == "Vol. 2 d"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 35, y = 0.1,
           label = "Study 2", hjust = 0, size = 10, color = "black") +
  annotate("text", x = Inf, y = Inf, label = "(d)", hjust = 1.1, vjust = 1.5, 
           size = 10, color = "black")

# See plot
plot.187

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB187VoluntSamplingRatesV5.png",
       plot = plot.187, width = 8, height = 10, dpi = 1300)

# Combine sampling rates --------------------------------------------------
# Combine the all data frames
combined_SR <- rbind(SR.V1.nd, SR.V1.d, SR.V2.d, SR.V3.1st.nd,
                     SR.V3.2nd.nd, SR_V3.2_nw, SR_V3.2_w)

# Capture the original order
congener_order <- combined_SR$congener %>% unique()

# Convert congener to factor with preserved order
combined_SR$congener <- factor(combined_SR$congener, levels = congener_order)

# Proceed with summarizing (do NOT use arrange here)
SR_averages_sd_cv <- combined_SR %>%
  group_by(congener) %>%
  summarise(
    n = sum(!is.na(`Sampling_Rate (m3/d)`)),
    Average_Sampling_Rate = if (n >= 3) mean(`Sampling_Rate (m3/d)`, na.rm = TRUE) else NA_real_,
    SD_Sampling_Rate = if (n >= 3) sd(`Sampling_Rate (m3/d)`, na.rm = TRUE) else NA_real_,
    CV_Sampling_Rate = if (n >= 3) (sd(`Sampling_Rate (m3/d)`, na.rm = TRUE) /
                                      mean(`Sampling_Rate (m3/d)`, na.rm = TRUE)) * 100 else NA_real_,
    Average_ko = if (n >= 3) mean(ko, na.rm = TRUE) else NA_real_
  ) %>%
  ungroup()

# Add ko values from averaging near values (4 above + center + 4 below)
# This is for future use of ko to determine Veff
SR_averages_sd_cv$Average_ko2 <- ifelse(
  is.na(SR_averages_sd_cv$Average_ko),
  zoo::rollapply(
    SR_averages_sd_cv$Average_ko,
    width = 9, # 4 above + center + 4 below
    FUN = function(x) mean(x, na.rm = TRUE),
    fill = NA,
    align = "center"
  ),
  SR_averages_sd_cv$Average_ko
)

# Add manually values to PCBs 207 to 209
SR_averages_sd_cv$Average_ko2[169] <- SR_averages_sd_cv$Average_ko2[167]
SR_averages_sd_cv$Average_ko2[170] <- SR_averages_sd_cv$Average_ko2[167]
SR_averages_sd_cv$Average_ko2[171] <- SR_averages_sd_cv$Average_ko2[167]

# Export results
write.csv(SR_averages_sd_cv,
          file = "Output/Data/csv/SamplingRates/Personal/PersonalAveSRV02.csv",
          row.names = FALSE)

# Plot the average and stdev
Plot.AV.SR <- ggplot(SR_averages_sd_cv, aes(x = congener, y = Average_Sampling_Rate)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_errorbar(aes(ymin = Average_Sampling_Rate,
                    ymax = Average_Sampling_Rate + SD_Sampling_Rate), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# See plot
print(Plot.AV.SR)

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/AvePersonalSRsV2.png",
       plot = Plot.AV.SR, width = 15, height = 5, dpi = 500)

# Plot the combined data with different colors for each group
Plot.SRs <- ggplot(combined_SR, aes(x = congener, y = `Sampling_Rate (m3/d)`,
                                    color = group)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "solid") +
  annotate("text", x = 2, y = 2, label = "Air Sampling Rate",
           hjust = -0.1, vjust = -1, color = "black", fontface = "bold") +  
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7)) +
  scale_color_discrete(name = "Participants")

# See plot
print(Plot.SRs)

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/SRsV02.png",
       plot = Plot.SRs, width = 15, height = 5, dpi = 500)

