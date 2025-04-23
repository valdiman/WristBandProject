## Script to calculate individual PCB "personal" sampling rates
# for silicone wristbands.
# nd non-dominant hand
# d dominant hand

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("RColorBrewer")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
}

# 3 volunteers, V1, V2, and V3
# Read data ---------------------------------------------------------------
data.V1 <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                             col_names = TRUE, col_types = NULL))
data.V2 <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.V3.1 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.V3.2 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet2",
                                  col_names = TRUE, col_types = NULL))
logKoa <- data.frame(read_excel("Data/logKoa.xlsx", sheet = "logKoa",
                                col_names = TRUE, col_types = NULL))
# ko from SamplingRates_ko.R file
ko <- read.csv("Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV1.csv")
# Select only ko [m/d]
ko <- ko[c(2,6)]

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
data.V1.pcbs <- data.V1[, c("time.day", "congeners", common_pcbs)]
data.V2.pcbs <- data.V2[, c("time.day", "congeners", common_pcbs)]
data.V3.1.pcbs <- data.V3.1[, c("time.day", "congeners", common_pcbs)]
data.V3.2.pcbs <- data.V3.2[, c("time.day", "congeners", common_pcbs)]

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
  Vwb <- 4.73 * 10^-6 # [m3]
  Awb <- 0.0054773 # [m2]
  # Calculate efective volume for static WBs
  veff_static.V1 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko.common$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.V1[1, 1]))
  # Compute concentration
  conc.V1 <- data.V1.2 / veff_static.V1
  # Calculate effective volume (Veff)
  subset_data <- data.V1.pcbs[4:13, 3:173]
  Veff.V1 <- t(apply(subset_data, 1, function(row) row / conc))
  # Add metadata to Veff.amanda and change format
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

SR.V1.nd <- data.frame(SR.V1.nd, group = "ParticipantV1.nd")
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
Awb.y = 0.0048707 # [m2] youth
SR.V1.nd$ko <- SR.V1.nd$`Sampling_Rate (m3/d)` / Awb.y # [m/d]

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

SR.V1.d <- data.frame(SR.V1.d, group = "ParticipantV1.d")
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
# Calculate ko from Amanda d
SR.V1.d$ko <- SR.V1.d$`Sampling_Rate (m3/d)` / Awb.s # [m/d]

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

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

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
  Vwb <- 4.73 * 10^-6 # [m3]
  Awb <- 0.0054773 # [m2]
  veff_stat.V2 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko.common$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.V2[1, 1]))
  # Compute concentration
  conc.V2 <- data.V2.2 / veff_stat.V2
  # Calculate effective volume (Veff)
  subset_data <- data.V2.pcbs[4:8, 3:173]
  Veff.V2 <- t(apply(subset_data, 1, function(row) row / conc.V2))
  # Add metadata to Veff.kay and change format
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

SR.V2.d <- data.frame(SR.V2.d, group = "ParticipantV2.d")
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
Awb.a = 0.0054773 # [m2] adult
SR.V2.d$ko <- SR.V2.d$`Sampling_Rate (m3/d)` / Awb.a # [m/d]

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

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

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
  Vwb <- 4.73e-6  # m3
  Awb <- 0.0054773  # m2
  # Prepare data
  data.V3.1.1 <- data.V3.1.pcbs[1:6, 3:173]
  logKwb_val <- 10^logKwb$logKwb
  ko_val <- ko.common$ko
  # Initialize lists to store results
  veff_stat_list <- list()
  conc_list <- list()
  Veff_list <- list()
  for (i in 1:6) {
    # Effective volume calculation
    veff_stat <- logKwb_val * Vwb * (1 - exp(-ko_val * Awb / Vwb / logKwb_val * data.V3.1[i, 1]))
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

SR.V3.1st.nd <- data.frame(SR.V3.1st.nd, group = "ParticipantV3.1st.nd")
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
Awb.a <- 0.0054773 # [m2]
SR.V3.1st.nd$ko <- SR.V3.1st.nd$`Sampling_Rate (m3/d)` / Awb.a # [m/d]

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

SR.V3.2nd.nd <- data.frame(SR.V3.2nd.nd, group = "ParticipantV3.2nd.nd")
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
Awb.a <- 0.0054773 # [m2]
SR.V3.2nd.nd$ko <- SR.V3.2nd.nd$`Sampling_Rate (m3/d)` / Awb.a # [m/d]

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

# V3 SR vs logKoa regression 1 ------------------------------------------
# (1) Average both d and nd
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
# Create a long dataframe combining SR.yau.1st. nd and SR.yau.2nd.nd
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
  annotate("text", x = min(sr.long.V3$logKoa) + 1.4, y = max(sr.long.V3$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.long.V3$logKoa) + 0.35, y = max(sr.long.V3$sr) * 1.13,
           label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
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
  annotate("text", x = 7.5, y = 7,
           label = paste("Vol. 3 (nd 1st week)"),size = 5) +
  annotate("text", x = 7.85, y = 6.7,
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

# Calculate personal sampling rate V3 2nd -------------------------------
# 3 WBs were not wiped and 3 WBs were wiped
# Dominant hand used here
# WBs were used to calculate PCB concentration
# Concentrations were calculated for each sampling day 
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Constants
  Vwb <- 4.73e-6  # m3
  Awb <- 0.0054773  # m2
  # Prepare data
  data.V3.2.1 <- data.V3.2.pcbs[1:3, 3:173]
  logKwb_val <- 10^logKwb$logKwb
  ko_val <- ko.common$ko
  # Initialize lists to store results
  veff_stat_list <- list()
  conc_list <- list()
  Veff_list <- list()
  for (i in 1:3) {
    # Effective volume calculation
    veff_stat <- logKwb_val * Vwb * (1 - exp(-ko_val * Awb / Vwb / logKwb_val * data.V3.2[i, 1]))
    veff_stat_list[[i]] <- veff_stat
    # Air concentration calculation
    conc <- data.V3.2.1[i, ] / veff_stat
    conc_list[[i]] <- conc
    # Effective volume (Veff)
    Veff <- data.V3.2.pcbs[i + 3, 3:173] / conc
    Veff_list[[i]] <- Veff
  }
  
  # Combine Veff results
  Veff.V3.2 <- do.call(rbind, Veff_list)
  Veff.V3.2 <- cbind(data.V3.2[4:9, 1:2], Veff.V3.2)
  Veff.V3.2 <- as.data.frame(Veff.V3.2)
  # Select w and nw
  Veff.V3.nw <- Veff.V3.2[1:3, 3:173]
  Veff.V3.nw.t <- Veff.V3.2[1:3, 1]
  Veff.V3.w <- Veff.V3.2[4:6, 3:173]
  Veff.V3.w.t <- Veff.V3.2[4:6, 1]
}

# Calculate sampling rate (SR) for nw and nw (m3/d)
# Create matrix for sampling rate (SR)
SR.V3.nw.d <- matrix(nrow = length(Veff.V3.nw[1,]), ncol = 3)

for(i in 1:length(SR.V3.nw.d[, 1])) {
  if (sum(!is.na(Veff.V3.nw[, i ]) & !is.infinite(Veff.V3.nw[, i])) == 3) {
    fit <- lm(Veff.V3.nw[, i] ~ 0 + Veff.V3.nw.t)
    SR.V3.nw.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.V3.nw.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.V3.nw.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.V3.nw.d[i, 1] <- 0
    SR.V3.nw.d[i, 2] <- 0
    SR.V3.nw.d[i, 3] <- 0
  }
}
  
SR.V3.nw.d <- data.frame(SR.V3.nw.d, group = "ParticipantV3.nw.d")
colnames(SR.V3.nw.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.V3.nw)[0, ])
SR.V3.nw.d <- cbind(congener, SR.V3.nw.d)

# Convert R2 and p-value to numeric
SR.V3.nw.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.V3.nw.d$`Sampling_Rate (m3/d)`)
SR.V3.nw.d$R2 <- as.numeric(SR.V3.nw.d$R2)
SR.V3.nw.d$p_value <- as.numeric(SR.V3.nw.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.V3.nw.d$R2 < 0.9 | SR.V3.nw.d$p_value > 0.05
SR.V3.nw.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.V3.nw.d$R2[mask] <- NA
SR.V3.nw.d$p_value[mask] <- NA
# Calculate ko from V3.2
Awb.a <- 0.0054773 # [m2]
SR.V3.nw.d$ko <- SR.V3.nw.d$`Sampling_Rate (m3/d)` / Awb.a # [m/d]

# Export results
write.csv(SR.V3.nw.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.V3.nw.d.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.V3.nw.d$congener <- factor(SR.V3.nw.d$congener,
                              levels = unique(SR.V3.nw.d$congener))

ggplot(SR.V3.nw.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# (2) Remove metadata from Veff.yau.w
Veff.yau.w.d.2 <- Veff.yau.w.d[, 3:175]

# Create matrix for sampling rate (SR)
SR.yau.w.d <- matrix(nrow = length(Veff.yau.w.d.2[1,]), ncol = 3)

for(i in 1:length(SR.yau.w.d[, 1])) {
  if (sum(!is.na(Veff.yau.w.d.2[, i ]) & !is.infinite(Veff.yau.w.d.2[, i])) == 3) {
    fit <- lm(Veff.yau.w.d.2[, i] ~ 0 + Veff.yau.w.d.t)
    SR.yau.w.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.w.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.w.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.w.d[i, 1] <- 0
    SR.yau.w.d[i, 2] <- 0
    SR.yau.w.d[i, 3] <- 0
  }
}

SR.yau.w.d <- data.frame(SR.yau.w.d, group = "ParticipantY.w.d")
colnames(SR.yau.w.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.yau.w.d.2)[0, ])
SR.yau.w.d <- cbind(congener, SR.yau.w.d)

# Convert R2 and p-value to numeric
SR.yau.w.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.yau.w.d$`Sampling_Rate (m3/d)`)
SR.yau.w.d$R2 <- as.numeric(SR.yau.w.d$R2)
SR.yau.w.d$p_value <- as.numeric(SR.yau.w.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.yau.w.d$R2 < 0.9 | SR.yau.w.d$p_value > 0.05
SR.yau.w.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.yau.w.d$R2[mask] <- NA
SR.yau.w.d$p_value[mask] <- NA

# Export results
write.csv(SR.yau.w.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.yau.w.d.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.yau.w.d$congener <- factor(SR.yau.w.d$congener,
                              levels = unique(SR.yau.w.d$congener))

ggplot(SR.yau.w.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Ya'u SR vs logKoa regression 2 ------------------------------------------
# Check difference btw w and nw
# (1) nw
sr.yau.nw <- as.data.frame(SR.yau.nw.d$`Sampling_Rate (m3/d)`, na.rm = TRUE)

sr.yau.nw$logkoa <- logKoa$logKoa
colnames(sr.yau.nw) <- c('sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.nw <- lm(log(sr.yau.nw$sr) ~ sr.yau.nw$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.nw)[1])  # exponentiate the intercept
b <- coef(model.yau.nw)[2]       # coefficient for logKoa
r2 <- summary(model.yau.nw)$r.squared

# plot
p.sr.yau.koa.nw <- ggplot(sr.yau.nw, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.3, y = 15,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 4) + 
  annotate("text", x = 6.75, y = 14.2,
           label = paste("R² = ", round(r2, 2)), size = 4) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.yau.koa.nw

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.nw.png",
       plot = p.sr.yau.koa.nw, width = 6, height = 6, dpi = 500)

# (2) w
sr.yau.w <- as.data.frame(SR.yau.w.d$`Sampling_Rate (m3/d)`, na.rm = TRUE)

sr.yau.w$logkoa <- logKoa$logKoa
colnames(sr.yau.w) <- c('sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.w <- lm(log(sr.yau.w$sr) ~ sr.yau.w$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.w)[1])  # exponentiate the intercept
b <- coef(model.yau.w)[2]       # coefficient for logKoa
r2 <- summary(model.yau.w)$r.squared

# plot
p.sr.yau.koa.w <- ggplot(sr.yau.w, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.3, y = 15,
           label = paste("sr = ", round(a, 3),
                         " * exp(", round(b, 2), " x log Koa)", sep = ""),
           size = 4) + 
  annotate("text", x = 6.75, y = 14.2,
           label = paste("R² = ", round(r2, 2)), size = 4) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.yau.koa.w

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.w.png",
       plot = p.sr.yau.koa.w, width = 6, height = 6, dpi = 500)

# Plot individual congeners -----------------------------------------------
# Combine plot
# Combine data, padding shorter dataset with NA
max_len <- max(length(Veff.amanda.d.t), length(Veff.amanda.nd.t),
               length(Veff.kay.d.t), length(Veff.yau.1st.nd.t),
               length(Veff.yau.2nd.nd.t), length(Veff.yau.w.d.t),
               length(Veff.yau.nw.d.t))

# PCBs 18+30
up.Amanda.d <- data.frame(
  time = c(Veff.amanda.d.t, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  veff = c(Veff.amanda.d$PCB18.30, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.Amanda.nd <- data.frame(
  time = c(Veff.amanda.nd.t, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  veff = c(Veff.amanda.nd$PCB18.30, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.kay.d <- data.frame(
  time = c(Veff.kay.d.t, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  veff = c(Veff.kay.d$PCB18.30, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  group = rep("Vol. 2 d")
)

up.yau.1st.nd <- data.frame(
  time = c(Veff.yau.1st.nd.t, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  veff = c(Veff.yau.1st.nd$PCB18.30, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  group = rep("Vol. 3 1st nd")
)

up.yau.2nd.nd <- data.frame(
  time = c(Veff.yau.2nd.nd.t, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  veff = c(Veff.yau.2nd.nd$PCB18.30, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  group = rep("Vol. 3 2nd nd")
)

up.yau.w.d <- data.frame(
  time = c(Veff.yau.w.d.t, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  veff = c(Veff.yau.w.d$PCB18.30, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  group = rep("Vol. 3 w d")
)

up.yau.nw.d <- data.frame(
  time = c(Veff.yau.nw.d.t, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  veff = c(Veff.yau.nw.d$PCB18.30, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  group = rep("Vol. 3 nw d")
)

combined_data <- rbind(up.Amanda.d, up.Amanda.nd, up.kay.d, up.yau.1st.nd,
                       up.yau.2nd.nd, up.yau.w.d, up.yau.nw.d)

combined_data <- rbind(up.Amanda.d, up.Amanda.nd, up.kay.d)

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
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
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
  annotate("text", x = 35, y = 0.05,
           label = "Study 2", hjust = 0, size = 10, color = "black")

# See plot
plot.18.30

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB18.30VoluntSamplingRates.png",
       plot = plot.18.30, width = 8, height = 10, dpi = 1300)

# PCB 52
up.Amanda.d <- data.frame(
  time = c(Veff.amanda.d.t, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  veff = c(Veff.amanda.d$PCB52, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.Amanda.nd <- data.frame(
  time = c(Veff.amanda.nd.t, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  veff = c(Veff.amanda.nd$PCB52, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.kay.d <- data.frame(
  time = c(Veff.kay.d.t, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  veff = c(Veff.kay.d$PCB52, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  group = rep("Vol. 2 d")
)

up.yau.1st.nd <- data.frame(
  time = c(Veff.yau.1st.nd.t, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  veff = c(Veff.yau.1st.nd$PCB52, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  group = rep("Vol. 3 1st nd")
)

up.yau.2nd.nd <- data.frame(
  time = c(Veff.yau.2nd.nd.t, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  veff = c(Veff.yau.2nd.nd$PCB52, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  group = rep("Vol. 3 2nd nd")
)

up.yau.w.d <- data.frame(
  time = c(Veff.yau.w.d.t, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  veff = c(Veff.yau.w.d$PCB52, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  group = rep("Vol. 3 w d")
)

up.yau.nw.d <- data.frame(
  time = c(Veff.yau.nw.d.t, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  veff = c(Veff.yau.nw.d$PCB52, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  group = rep("Vol. 3 nw d")
)

combined_data <- rbind(up.Amanda.d, up.Amanda.nd, up.kay.d)

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
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
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
           label = "Study 2", hjust = 0, size = 10, color = "black")

# See plot
plot.52

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB52VoluntSamplingRates.png",
       plot = plot.52, width = 8, height = 10, dpi = 1300)

# PCB 118
up.Amanda.d <- data.frame(
  time = c(Veff.amanda.d.t, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  veff = c(Veff.amanda.d$PCB118, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.Amanda.nd <- data.frame(
  time = c(Veff.amanda.nd.t, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  veff = c(Veff.amanda.nd$PCB118, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.kay.d <- data.frame(
  time = c(Veff.kay.d.t, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  veff = c(Veff.kay.d$PCB118, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  group = rep("Vol. 2 d")
)

up.yau.1st.nd <- data.frame(
  time = c(Veff.yau.1st.nd.t, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  veff = c(Veff.yau.1st.nd$PCB118, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  group = rep("Vol. 3 1st nd")
)

up.yau.2nd.nd <- data.frame(
  time = c(Veff.yau.2nd.nd.t, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  veff = c(Veff.yau.2nd.nd$PCB118, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  group = rep("Vol. 3 2nd nd")
)

up.yau.w.d <- data.frame(
  time = c(Veff.yau.w.d.t, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  veff = c(Veff.yau.w.d$PCB118, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  group = rep("Vol. 3 w d")
)

up.yau.nw.d <- data.frame(
  time = c(Veff.yau.nw.d.t, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  veff = c(Veff.yau.nw.d$PCB118, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  group = rep("Vol. 3 nw d")
)

combined_data <- rbind(up.Amanda.d, up.Amanda.nd, up.kay.d)

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
ggsave("Output/Plots/SamplingRates/Personal/PCB118VoluntSamplingRates.png",
       plot = plot.118, width = 8, height = 10, dpi = 1300)

# PCB 187
up.Amanda.d <- data.frame(
  time = c(Veff.amanda.d.t, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  veff = c(Veff.amanda.d$PCB187, rep(NA, max_len - length(Veff.amanda.d[, 1]))),
  group = rep("Vol. 1 d")
)

up.Amanda.nd <- data.frame(
  time = c(Veff.amanda.nd.t, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  veff = c(Veff.amanda.nd$PCB187, rep(NA, max_len - length(Veff.amanda.nd[, 1]))),
  group = rep("Vol. 1 nd")
)

up.kay.d <- data.frame(
  time = c(Veff.kay.d.t, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  veff = c(Veff.kay.d$PCB187, rep(NA, max_len - length(Veff.kay.d[, 1]))),
  group = rep("Vol. 2 d")
)

up.yau.1st.nd <- data.frame(
  time = c(Veff.yau.1st.nd.t, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  veff = c(Veff.yau.1st.nd$PCB187, rep(NA, max_len - length(Veff.yau.1st.nd[, 1]))),
  group = rep("Vol. 3 1st nd")
)

up.yau.2nd.nd <- data.frame(
  time = c(Veff.yau.2nd.nd.t, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  veff = c(Veff.yau.2nd.nd$PCB187, rep(NA, max_len - length(Veff.yau.2nd.nd[, 1]))),
  group = rep("Vol. 3 2nd nd")
)

up.yau.w.d <- data.frame(
  time = c(Veff.yau.w.d.t, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  veff = c(Veff.yau.w.d$PCB187, rep(NA, max_len - length(Veff.yau.w.d[, 1]))),
  group = rep("Vol. 3 w d")
)

up.yau.nw.d <- data.frame(
  time = c(Veff.yau.nw.d.t, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  veff = c(Veff.yau.nw.d$PCB187, rep(NA, max_len - length(Veff.yau.nw.d[, 1]))),
  group = rep("Vol. 3 nw d")
)

combined_data <- rbind(up.Amanda.d, up.Amanda.nd, up.kay.d)

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
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
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
  annotate("text", x = 35, y = 0.6,
           label = "Study 2", hjust = 0, size = 10, color = "black")

# See plot
plot.187

# Save plot
ggsave("Output/Plots/SamplingRates/Personal/PCB187VoluntSamplingRates.png",
       plot = plot.187, width = 8, height = 10, dpi = 1300)

# Combine sampling rates --------------------------------------------------
# Combine the all data frames
combined_SR <- rbind(SR.amanda.nd, SR.amanda.d, SR.kay.d, SR.yau.1st.nd,
                     SR.yau.2nd.nd, SR.yau.nw.d, SR.yau.w.d)

# Look at SR and variability
SR_averages_sd_cv <- combined_SR %>%
  group_by(congener) %>%
  summarise(
    Average_Sampling_Rate = mean(`Sampling_Rate (m3/d)`, na.rm = TRUE),
    SD_Sampling_Rate = sd(`Sampling_Rate (m3/d)`, na.rm = TRUE),
    CV_Sampling_Rate = (sd(`Sampling_Rate (m3/d)`,
                           na.rm = TRUE) / mean(`Sampling_Rate (m3/d)`, na.rm = TRUE)) * 100
  )

# Remove PCB (rows) with only one SR measurement, i.e., SD and CV = NA
SR_averages_sd_cv <- SR_averages_sd_cv %>%
  filter(!is.na(SD_Sampling_Rate))

# Export results
write.csv(SR_averages_sd_cv,
          file = "Output/Data/csv/SamplingRates/Personal/PersonalAveSRV01.csv",
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
ggsave("Output/Plots/SamplingRates/Personal/AvePersonalSRs.png",
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
ggsave("Output/Plots/SamplingRates/Personal/SRsV01.png",
       plot = Plot.SRs, width = 15, height = 5, dpi = 500)

# PCB profiles ------------------------------------------------------------
# Profiles should be created using concentration not mass
# This is just to review the PCB distributions

# Amanda ------------------------------------------------------------------
# Select data 
stat.amanda <- t(data.frame(data.amanda.2))
worn.amanda.R <- data.amanda[8, 3:175]
worn.amanda.L <- data.amanda[13, 3:175]
value.amanda <- rbind(stat.amanda, worn.amanda.R, worn.amanda.L)
tmp <- rowSums(value.amanda, na.rm = TRUE)
profile.amanda <- sweep(value.amanda, 1, tmp, FUN = "/")
profile.amanda <- t(profile.amanda)
profile_matrix <- as.matrix(profile.amanda)
profile.amanda <- data.frame(RowNames = rownames(profile_matrix),
                             profile_matrix)
rownames(profile.amanda) <- NULL
colnames(profile.amanda) <- c("congeners", "Air", "ParticipantA.r",
                              "ParticipantA.l")
profile.amanda$congeners <- factor(profile.amanda$congeners,
                                   levels = unique(profile.amanda$congeners))

# Reshape the data frame to long format
profile_long <- profile.amanda %>%
  pivot_longer(cols = c(Air, ParticipantA.r, ParticipantA.l),
               names_to = "Participant",
               values_to = "Value")

# Define color
palette <- brewer.pal(3, "Set1")

# Profile plot
Plot.prof.amanda <- ggplot(profile_long, aes(x = congeners, y = Value,
                                             fill = Participant)) +
  geom_bar(stat = "identity", position = "dodge", width = 1, alpha = 0.8) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60,
                                   hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  scale_fill_manual(name = "Samples", values = palette)

# see plot
print(Plot.prof.amanda)

# 1:1 plots
threshold <- 0.005  # Define the threshold for labeling

i <- ggplot(profile.amanda, aes(x = ParticipantA.l,
                                y = ParticipantA.r, label = congeners)) +
  geom_point() +
  geom_text(data = subset(profile.amanda,
                          abs(ParticipantA.l - ParticipantA.r) > threshold),
            size = 3, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("WB Personal A (l)"))) +
  xlab(expression(bold("WB Personal A (r)")))

ii <- ggplot(profile.amanda, aes(x = Air, y = ParticipantA.l,
                                 label = congeners)) +
  geom_point() +
  geom_text(data = subset(profile.amanda,
                          abs(Air - ParticipantA.l) > threshold),
            size = 3, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air"))) +
  xlab(expression(bold("WB Personal A (l)")))

iii <- ggplot(profile.amanda, aes(x = Air, y = ParticipantA.r,
                                  label = congeners)) +
  geom_point() +
  geom_text(data = subset(profile.amanda,
                          abs(Air - ParticipantA.r) > threshold),
            size = 3, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air"))) +
  xlab(expression(bold("WB Personal A (r)")))

combined_plot <- grid.arrange(i, ii, iii, nrow = 1)

# Kay ---------------------------------------------------------------------
# Select data 
stat.kay <- t(data.frame(data.kay.2))
worn.kay.R <- data.kay[8, 3:175]
value.kay <- rbind(stat.kay, worn.kay.R)
tmp <- rowSums(value.kay, na.rm = TRUE)
profile.kay <- sweep(value.kay, 1, tmp, FUN = "/")
profile.kay <- t(profile.kay)
profile_matrix <- as.matrix(profile.kay)
profile.kay <- data.frame(RowNames = rownames(profile_matrix),
                             profile_matrix)
rownames(profile.kay) <- NULL
colnames(profile.kay) <- c("congeners", "Air", "ParticipantK.r")
profile.kay$congeners <- factor(profile.kay$congeners,
                                   levels = unique(profile.kay$congeners))

# Reshape the data frame to long format
profile_long <- profile.kay %>%
  pivot_longer(cols = c(Air, ParticipantK.r),
               names_to = "Samples",
               values_to = "Value")

# Profile plot
Plot.prof.kay <- ggplot(profile_long, aes(x = congeners, y = Value,
                                             fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 1, alpha = 0.8) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60,
                                   hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB")))

# see plot
print(Plot.prof.kay)

# 1:1 plots
threshold <- 0.002  # Define the threshold for labeling

i <- ggplot(profile.kay, aes(x = Air,
                                y = ParticipantK.r, label = congeners)) +
  geom_point() +
  geom_text(data = subset(profile.kay,
                          abs(Air - ParticipantK.r) > threshold),
            size = 3, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air"))) +
  xlab(expression(bold("WB Personal K")))

# see plot
print(i)


# Ya'u --------------------------------------------------------------------
tmp <- rowSums(data.yau[, 4:176], na.rm = TRUE)
profile.yau <- sweep(data.yau[, 4:176], 1, tmp, FUN = "/")
profile.yau <- cbind(data.yau$congeners, profile.yau)
profile.yau <- cbind(data.yau$week, profile.yau)
names(profile.yau)[1] <- "week"
# Select week and days
# Week 1
profile.yau.1st <- profile.yau[profile.yau$week == 1, ]
profile.yau.1st <- t(profile.yau.1st)
profile.yau.1st <- profile.yau.1st[-1, ]
colnames(profile.yau.1st) <- profile.yau.1st[1,]
profile.yau.1st <- profile.yau.1st[-1, ]
rownames_data <- rownames(profile.yau.1st)
rownames(profile.yau.1st) <- NULL
profile.yau.1st <- cbind(Row_Name = rownames_data, profile.yau.1st)
colnames(profile.yau.1st)[1] <- "congeners"
colnames(profile.yau.1st)[2] <- "Air.day1"
colnames(profile.yau.1st)[3] <- "Air.day3"
colnames(profile.yau.1st)[4] <- "Air.day5"
colnames(profile.yau.1st)[5] <- "WB Personal Y (1st) day 1"
colnames(profile.yau.1st)[6] <- "WB Personal Y (1st) day 3"
colnames(profile.yau.1st)[7] <- "WB Personal Y (1st) day 5"

profile_long <- profile.yau.1st %>%
  as.data.frame() %>%
  pivot_longer(cols = -congeners, 
               names_to = "Samples", 
               values_to = "Value")

# Day 1
selected_samples <- c("Air.day1", "WB Personal Y (1st) day 1")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, max(profile_long_filtered$Value, na.rm = TRUE) * 1.1) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
prof.yau.1st <- data.frame(profile.yau.1st)
prof.yau.1st[, 2:7] <- lapply(prof.yau.1st[, 2:7],
                              as.numeric)
threshold <- 0.005  # Define the threshold for labeling

Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day1,
                                   y = WB.Personal.Y..1st..day.1,
                                   label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day1 - WB.Personal.Y..1st..day.1) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 1"))) +
  xlab(expression(bold("WB Personal Y (1st) day 1")))

# see plot
print(Plot.prof.yau)

# Day 3
selected_samples <- c("Air.day3", "WB Personal Y (1st) day 3")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day3,
                                          y = WB.Personal.Y..1st..day.3,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day3 - WB.Personal.Y..1st..day.3) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 3"))) +
  xlab(expression(bold("WB Personal Y (1st) Day 3")))

# see plot
print(Plot.prof.yau)

# Day 5
selected_samples <- c("Air.day5", "WB Personal Y (1st) day 5")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day5,
                                          y = WB.Personal.Y..1st..day.5,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day5 - WB.Personal.Y..1st..day.5) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air day 5"))) +
  xlab(expression(bold("WB Personal Y (1st) day 5")))

# see plot
print(Plot.prof.yau)

# Week 2
profile.yau.2nd <- profile.yau[profile.yau$week == 2, ]
profile.yau.2nd <- t(profile.yau.2nd)
profile.yau.2nd <- profile.yau.2nd[-1, ]
colnames(profile.yau.2nd) <- profile.yau.2nd[1,]
profile.yau.2nd <- profile.yau.2nd[-1, ]
rownames_data <- rownames(profile.yau.2nd)
rownames(profile.yau.2nd) <- NULL
profile.yau.2nd <- cbind(Row_Name = rownames_data, profile.yau.2nd)
colnames(profile.yau.2nd)[1] <- "congeners"
colnames(profile.yau.2nd)[2] <- "Air.day1"
colnames(profile.yau.2nd)[3] <- "Air.day3"
colnames(profile.yau.2nd)[4] <- "Air.day5"
colnames(profile.yau.2nd)[5] <- "WB Personal Y (2nd) day 1"
colnames(profile.yau.2nd)[6] <- "WB Personal Y (2nd) day 3"
colnames(profile.yau.2nd)[7] <- "WB Personal Y (2nd) day 5"

profile_long <- profile.yau.2nd %>%
  as.data.frame() %>%
  pivot_longer(cols = -congeners, 
               names_to = "Samples", 
               values_to = "Value")

# Day 1
selected_samples <- c("Air.day1", "WB Personal Y (2nd) day 1")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
prof.yau.2nd <- data.frame(profile.yau.2nd)
prof.yau.2nd[, 2:7] <- lapply(prof.yau.2nd[, 2:7],
                              as.numeric)
threshold <- 0.005  # Define the threshold for labeling

Plot.prof.yau <- ggplot(prof.yau.2nd, aes(x = Air.day1,
                                          y = WB.Personal.Y..2nd..day.1,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.2nd,
                          abs(Air.day1 - WB.Personal.Y..2nd..day.1) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 1"))) +
  xlab(expression(bold("WB Personal Y (2nd) day 1")))

# see plot
print(Plot.prof.yau)

# Day 3
selected_samples <- c("Air.day3", "WB Personal Y (2nd) day 3")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.2nd, aes(x = Air.day3,
                                          y = WB.Personal.Y..2nd..day.3,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.2nd,
                          abs(Air.day3 - WB.Personal.Y..2nd..day.3) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 3"))) +
  xlab(expression(bold("WB Personal Y (2nd) Day 3")))

# see plot
print(Plot.prof.yau)

# Day 5
selected_samples <- c("Air.day5", "WB Personal Y (2nd) day 5")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.2nd, aes(x = Air.day5,
                                          y = WB.Personal.Y..2nd..day.5,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.2nd,
                          abs(Air.day5 - WB.Personal.Y..2nd..day.5) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 5"))) +
  xlab(expression(bold("WB Personal Y (2nd) Day 5")))

# see plot
print(Plot.prof.yau)

# Ya'u nw & w
tmp <- rowSums(data.yau2[, 5:177], na.rm = TRUE)
profile.yau <- sweep(data.yau2[, 5:177], 1, tmp, FUN = "/")
profile.yau <- cbind(data.yau2$congeners, profile.yau)
profile.yau <- cbind(data.yau2$wiped, profile.yau)
profile.yau <- cbind(data.yau2$static, profile.yau)
names(profile.yau)[1] <- "static"
names(profile.yau)[2] <- "wiped"
names(profile.yau)[3] <- "sample"
# Select static and days
# profile.yau.st <- profile.yau[profile.yau$static == 'y', ]
profile.yau <- t(profile.yau)
profile.yau <- profile.yau[-c(1:3), ]
rownames_data <- rownames(profile.yau)
rownames(profile.yau) <- NULL
profile.yau <- cbind(Row_Name = rownames_data, profile.yau)
{
  colnames(profile.yau)[1] <- "congeners"
  colnames(profile.yau)[2] <- "Air.day1"
  colnames(profile.yau)[3] <- "Air.day3"
  colnames(profile.yau)[4] <- "Air.day5"
  colnames(profile.yau)[5] <- "WB Personal Y (nw) day 1"
  colnames(profile.yau)[6] <- "WB Personal Y (nw) day 3"
  colnames(profile.yau)[7] <- "WB Personal Y (nw) day 5"
  colnames(profile.yau)[8] <- "WB Personal Y (w) day 1"
  colnames(profile.yau)[9] <- "WB Personal Y (w) day 3"
  colnames(profile.yau)[10] <- "WB Personal Y (w) day 5"
}

profile_long <- profile.yau %>%
  as.data.frame() %>%
  pivot_longer(cols = -congeners, 
               names_to = "Samples", 
               values_to = "Value")

# Day 1
selected_samples <- c("Air.day1", "WB Personal Y (nw) day 1",
                      "WB Personal Y (w) day 1")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, max(profile_long_filtered$Value, na.rm = TRUE) * 1.1) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
prof.yau.1st <- data.frame(profile.yau.1st)
prof.yau.1st[, 2:7] <- lapply(prof.yau.1st[, 2:7],
                              as.numeric)
threshold <- 0.005  # Define the threshold for labeling

Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day1,
                                          y = WB.Personal.Y..1st..day.1,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day1 - WB.Personal.Y..1st..day.1) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 1"))) +
  xlab(expression(bold("WB Personal Y (1st) day 1")))

# see plot
print(Plot.prof.yau)

# Day 3
selected_samples <- c("Air.day3", "WB Personal Y (nw) day 3",
                      "WB Personal Y (w) day 3")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day3,
                                          y = WB.Personal.Y..1st..day.3,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day3 - WB.Personal.Y..1st..day.3) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air Day 3"))) +
  xlab(expression(bold("WB Personal Y (1st) Day 3")))

# see plot
print(Plot.prof.yau)

# Day 5
selected_samples <- c("Air.day5", "WB Personal Y (nw) day 5",
                      "WB Personal Y (w) day 5")
profile_long_filtered <- profile_long %>%
  filter(Samples %in% selected_samples)
profile_long_filtered$Value <- as.numeric(profile_long_filtered$Value)
profile_long_filtered$congeners <- factor(profile_long_filtered$congeners,
                                          levels = unique(profile_long_filtered$congeners))

# Bar plot
Plot.prof.yau <- ggplot(profile_long_filtered,
                        aes(x = congeners, y = Value, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12,
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13))

# see plot
print(Plot.prof.yau)

# 1:1 plot
Plot.prof.yau <- ggplot(prof.yau.1st, aes(x = Air.day5,
                                          y = WB.Personal.Y..1st..day.5,
                                          label = congeners)) +
  geom_point(size = 3) +
  geom_text(data = subset(prof.yau.1st,
                          abs(Air.day5 - WB.Personal.Y..1st..day.5) > threshold),
            size = 5, vjust = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0, 0.15) +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 10/10,
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  ylab(expression(bold("Air day 5"))) +
  xlab(expression(bold("WB Personal Y (1st) day 5")))

# see plot
print(Plot.prof.yau)

