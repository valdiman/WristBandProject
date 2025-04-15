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

# Read data from excel ----------------------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                             col_names = TRUE, col_types = NULL))
data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau2 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet2",
                                  col_names = TRUE, col_types = NULL))
logKoa <- data.frame(read_excel("Data/logKoa.xlsx", sheet = "logKoa",
                                col_names = TRUE, col_types = NULL))

# Calculate personal sampling rate Amanda ---------------------------------
# WBs were used to calculate PCB concentration
# Both hands (d and nd)
# triplicates for 4.3 days were deployed
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.amanda.1 <- data.amanda[1:3, ]
  # Average 3 WBs. NA values not included in the calculations
  data.amanda.2 <- colMeans(data.amanda.1[, 3:175], na.rm = TRUE)
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc <- data.amanda.2/(0.5*data.amanda[1,1])
  # Calculate effective volume (Veff)
  subset_data <- data.amanda[4:13, 3:175]
  Veff.amanda <- t(apply(subset_data, 1, function(row) row / conc))
  # Add metadata to Veff.amanda and change format
  Veff.amanda <- cbind(data.amanda[4:13, 2], data.amanda[4:13, 1], Veff.amanda)
  # Transform to data.frame
  Veff.amanda <- as.data.frame(Veff.amanda)
  # Add names to first 2 columns
  colnames(Veff.amanda)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.amanda[, 2:175] <- apply(Veff.amanda[, 2:175], 2, as.numeric)
  # Select right, remove metadata
  Veff.amanda.nd <- Veff.amanda[1:5, 3:175]
  # Select time
  Veff.amanda.nd.t <- Veff.amanda[1:5, 2]
  # Select left, remove metadata
  Veff.amanda.d <- Veff.amanda[6:10, 3:175]
  # Select time
  Veff.amanda.d.t <- Veff.amanda[6:10, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.amanda.nd <- matrix(nrow = length(Veff.amanda.nd[1,]), ncol = 3)

for(i in 1:length(SR.amanda.nd[, 1])) {
  if (length(unique(Veff.amanda.nd[, i])) >= 3) {
    fit <- lm(Veff.amanda.nd[, i] ~ 0 + Veff.amanda.nd.t)
    SR.amanda.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.amanda.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.amanda.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.amanda.nd[i, 1] <- 0
    SR.amanda.nd[i, 2] <- 0
    SR.amanda.nd[i, 3] <- 0
  }
}

SR.amanda.nd <- data.frame(SR.amanda.nd, group = "ParticipantA.nd")
colnames(SR.amanda.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.amanda.nd)[0, ])
SR.amanda.nd <- cbind(congener, SR.amanda.nd)

# Convert R2 and p-value to numeric
SR.amanda.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.amanda.nd$`Sampling_Rate (m3/d)`)
SR.amanda.nd$R2 <- as.numeric(SR.amanda.nd$R2)
SR.amanda.nd$p_value <- as.numeric(SR.amanda.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.amanda.nd$R2 < 0.9 | SR.amanda.nd$p_value > 0.05
SR.amanda.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.amanda.nd$R2[mask] <- NA
SR.amanda.nd$p_value[mask] <- NA

# Export results
write.csv(SR.amanda.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.amanda.nd.csv", row.names = FALSE)

# Create matrix for sampling rate (SR)
SR.amanda.d <- matrix(nrow = length(Veff.amanda.d[1,]), ncol = 3)

for(i in 1:length(SR.amanda.d[, 1])) {
  if (length(unique(Veff.amanda.d[, i])) >= 3) {
    fit <- lm(Veff.amanda.d[, i] ~ 0 + Veff.amanda.d.t)
    SR.amanda.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.amanda.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.amanda.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.amanda.d[i, 1] <- 0
    SR.amanda.d[i, 2] <- 0
    SR.amanda.d[i, 3] <- 0
  }
}

SR.amanda.d <- data.frame(SR.amanda.d, group = "ParticipantA.d")
colnames(SR.amanda.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.amanda.d)[0, ])
SR.amanda.d <- cbind(congener, SR.amanda.d)

# Convert R2 and p-value to numeric
SR.amanda.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.amanda.d$`Sampling_Rate (m3/d)`)
SR.amanda.d$R2 <- as.numeric(SR.amanda.d$R2)
SR.amanda.d$p_value <- as.numeric(SR.amanda.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.amanda.d$R2 < 0.9 | SR.amanda.d$p_value > 0.05
SR.amanda.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.amanda.d$R2[mask] <- NA
SR.amanda.d$p_value[mask] <- NA

# Export results
write.csv(SR.amanda.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.amanda.d.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.amanda.nd$congener <- factor(SR.amanda.nd$congener,
                            levels = unique(SR.amanda.nd$congener))
# Plot with legend
ggplot(SR.amanda.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
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
SR.amanda.d$congener <- factor(SR.amanda.d$congener,
                               levels = unique(SR.amanda.d$congener))
# Plot with legend
ggplot(SR.amanda.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Amanda SR vs logKoa regression ------------------------------------------
# (1) Average both d and nd
sr.ave.amanda <- as.data.frame(rowMeans(cbind(SR.amanda.d$`Sampling_Rate (m3/d)`, 
                                SR.amanda.nd$`Sampling_Rate (m3/d)`), 
                          na.rm = TRUE))

sr.ave.amanda$logkoa <- logKoa$logKoa
colnames(sr.ave.amanda) <- c('ave_sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.amanda.1 <- lm(log(sr.ave.amanda$ave_sr) ~ sr.ave.amanda$logKoa)

# Get the coefficients
a <- exp(coef(model.amanda.1)[1])  # exponentiate the intercept
b <- coef(model.amanda.1)[2]       # coefficient for logKoa
r2 <- summary(model.amanda.1)$r.squared

# plot
p.sr.amanda.koa.1 <- ggplot(sr.ave.amanda, aes(x = logKoa, y = ave_sr)) +
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

p.sr.amanda.koa.1

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Amanda_logKoa1.png", plot = p.sr.amanda.koa.1,
       width = 6, height = 6, dpi = 500)

# (2) Individual values
# Create a long dataframe combining SR.amanda.d and SR.amanda.nd
sr.long.amanda <- data.frame(
  sr = c(SR.amanda.d$`Sampling_Rate (m3/d)`, SR.amanda.nd$`Sampling_Rate (m3/d)`),
  logKoa = rep(logKoa$logKoa, 2)  # Repeat logKoa values for both d and nd
)

# Remove any NA values
sr.long.amanda <- na.omit(sr.long.amanda)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.amanda.2 <- lm(log(sr.long.amanda$sr) ~ sr.long.amanda$logKoa)

# Get the coefficients
a <- exp(coef(model.amanda.2)[1])  # exponentiate the intercept
b <- coef(model.amanda.2)[2]       # coefficient for logKoa
r2 <- summary(model.amanda.2)$r.squared

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

# Plot
p.sr.amanda.koa.2 <- ggplot(sr.long.amanda, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.long.amanda$logKoa) + 0.6, y = max(sr.long.amanda$sr) * 1.2,
           label = paste("Vol. 1 (d & nd)"),size = 5) +
  annotate("text", x = min(sr.long.amanda$logKoa) + 1.5, y = max(sr.long.amanda$sr) * 1.15,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.long.amanda$logKoa) + 0.35, y = max(sr.long.amanda$sr) * 1.1,
           label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12))

p.sr.amanda.koa.2

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Amanda_logKoa2.png", plot = p.sr.amanda.koa.2,
       width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate Kay ------------------------------------
# WBs were used to calculate PCB concentration
# Dominant hand (d)
# triplicates for 4.27 days were deployed
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.kay.1 <- data.kay[1:3,]
  # Average 3 WBs
  data.kay.2 <- colMeans(data.kay.1[, 3:175])
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc <- data.kay.2/(0.5*data.kay[1,1])
  # Calculate effective volume (Veff)
  subset_data <- data.kay[4:8, 3:175]
  Veff.kay <- t(apply(subset_data, 1, function(row) row / conc))
  # Add metadata to Veff.kay and change format
  Veff.kay <- cbind(data.kay[4:8, 2], data.kay[4:8, 1], Veff.kay)
  # Transform to data.frame
  Veff.kay <- as.data.frame(Veff.kay)
  # Add names to first 2 columns
  colnames(Veff.kay)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.kay[, 2:175] <- apply(Veff.kay[, 2:175], 2, as.numeric)
  # Select right, remove metadata
  Veff.kay.d <- Veff.kay[1:5, 3:175]
  # Select time
  Veff.kay.d.t <- Veff.kay[1:5, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.kay.d <- matrix(nrow = length(Veff.kay.d[1,]), ncol = 3)

for(i in 1:length(SR.kay.d[, 1])) {
  if (length(unique(Veff.kay.d[, i])) >= 3) {
    fit <- lm(Veff.kay.d[, i] ~ 0 + Veff.kay.d.t)
    SR.kay.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.kay.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.kay.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.kay.d[i, 1] <- 0
    SR.kay.d[i, 2] <- 0
    SR.kay.d[i, 3] <- 0
  }
}

SR.kay.d <- data.frame(SR.kay.d, group = "ParticipantK.d")
colnames(SR.kay.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.kay.d)[0, ])
SR.kay.d <- cbind(congener, SR.kay.d)

# Convert R2 and p-value to numeric
SR.kay.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.kay.d$`Sampling_Rate (m3/d)`)
SR.kay.d$R2 <- as.numeric(SR.kay.d$R2)
SR.kay.d$p_value <- as.numeric(SR.kay.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.kay.d$R2 < 0.9 | SR.kay.d$p_value > 0.05
SR.kay.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.kay.d$R2[mask] <- NA
SR.kay.d$p_value[mask] <- NA

# Export results
write.csv(SR.kay.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.kay.d.csv", row.names = FALSE)

# Plot
# Organize PCB names
SR.kay.d$congener <- factor(SR.kay.d$congener,
                               levels = unique(SR.kay.d$congener))

ggplot(SR.kay.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Kay SR vs logKoa regression ---------------------------------------------
# Create a long dataframe combining SR.kay.d
sr.kay <- data.frame(
  sr = SR.kay.d$`Sampling_Rate (m3/d)`,
  logKoa = logKoa$logKoa)

# Remove any NA values
sr.kay <- na.omit(sr.kay)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.kay <- lm(log(sr.kay$sr) ~ sr.kay$logKoa)

# Get the coefficients
a <- exp(coef(model.kay)[1])  # exponentiate the intercept
b <- coef(model.kay)[2]       # coefficient for logKoa
r2 <- summary(model.kay)$r.squared

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

# Plot
p.sr.kay.koa <- ggplot(sr.kay, aes(x = logKoa, y = sr)) +
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

p.sr.kay.koa

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Kay_logKoa.png", plot = p.sr.kay.koa,
       width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate Ya'u -----------------------------------
# WBs were used to calculate PCB concentration
# Non-dominant hand (nd)
# Concentrations were calculated for each sampling day 
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.yau.1 <- data.yau[1:6, 4:176]
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  time <- data.yau[1:6, 1]
  conc <- t(sweep(data.yau.1, 1, 0.5 * time, "/"))
  # get WB mass
  mass.WD <- data.yau[7:12, 4:176]
  Veff.yau <- mass.WD/t(conc)
  # Add metadata to Veff.amanda and change format
  Veff.yau <- cbind(data.yau[7:12, 3], data.yau[7:12, 1], Veff.yau)
  # Transform to data.frame
  Veff.yau <- as.data.frame(Veff.yau)
  # Add names to first 2 columns
  colnames(Veff.yau)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.yau[, 2:175] <- apply(Veff.yau[, 2:175], 2, as.numeric)
  # Select 1st week, remove metadata
  Veff.yau.1st.nd <- Veff.yau[1:3, 3:175]
  # Select time
  Veff.yau.1st.nd.t <- Veff.yau[1:3, 2]
  # Select left, remove metadata
  Veff.yau.2nd.nd <- Veff.yau[4:6, 3:175]
  # Select time
  Veff.yau.2nd.nd.t <- Veff.yau[4:6, 2]
}

# Calculate sampling rate (SR) for 1st and 2nd weeks (m3/d)
# Create matrix for sampling rate (SR)
SR.yau.1st.nd <- matrix(nrow = length(Veff.yau.1st.nd[1,]), ncol = 3)

for(i in 1:length(SR.yau.1st.nd[, 1])) {
  if (sum(!is.na(Veff.yau.1st.nd[, i ]) & !is.infinite(Veff.yau.1st.nd[, i])) == 3) {
    fit <- lm(Veff.yau.1st.nd[, i] ~ 0 + Veff.yau.1st.nd.t)
    SR.yau.1st.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.1st.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.1st.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.1st.nd[i, 1] <- 0
    SR.yau.1st.nd[i, 2] <- 0
    SR.yau.1st.nd[i, 3] <- 0
  }
}

SR.yau.1st.nd <- data.frame(SR.yau.1st.nd, group = "ParticipantY.1st.nd")
colnames(SR.yau.1st.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.yau.1st.nd)[0, ])
SR.yau.1st.nd <- cbind(congener, SR.yau.1st.nd)

# Convert R2 and p-value to numeric
SR.yau.1st.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.yau.1st.nd$`Sampling_Rate (m3/d)`)
SR.yau.1st.nd$R2 <- as.numeric(SR.yau.1st.nd$R2)
SR.yau.1st.nd$p_value <- as.numeric(SR.yau.1st.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.yau.1st.nd$R2 < 0.9 | SR.yau.1st.nd$p_value > 0.05
SR.yau.1st.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.yau.1st.nd$R2[mask] <- NA
SR.yau.1st.nd$p_value[mask] <- NA

# Export results
write.csv(SR.yau.1st.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.yau.1st.nd.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.yau.1st.nd$congener <- factor(SR.yau.1st.nd$congener,
                            levels = unique(SR.yau.1st.nd$congener))

ggplot(SR.yau.1st.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
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
SR.yau.2nd.nd <- matrix(nrow = length(Veff.yau.2nd.nd[1,]), ncol = 3)

for(i in 1:length(SR.yau.2nd.nd[, 1])) {
  if (sum(!is.na(Veff.yau.2nd.nd[, i ]) & !is.infinite(Veff.yau.2nd.nd[, i])) == 3) {
    fit <- lm(Veff.yau.2nd.nd[, i] ~ 0 + Veff.yau.2nd.nd.t)
    SR.yau.2nd.nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.2nd.nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.2nd.nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.2nd.nd[i, 1] <- 0
    SR.yau.2nd.nd[i, 2] <- 0
    SR.yau.2nd.nd[i, 3] <- 0
  }
}

SR.yau.2nd.nd <- data.frame(SR.yau.2nd.nd, group = "ParticipantY.2nd.nd")
colnames(SR.yau.2nd.nd) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.yau.2nd.nd)[0, ])
SR.yau.2nd.nd <- cbind(congener, SR.yau.2nd.nd)

# Convert R2 and p-value to numeric
SR.yau.2nd.nd$`Sampling_Rate (m3/d)` <- as.numeric(SR.yau.2nd.nd$`Sampling_Rate (m3/d)`)
SR.yau.2nd.nd$R2 <- as.numeric(SR.yau.2nd.nd$R2)
SR.yau.2nd.nd$p_value <- as.numeric(SR.yau.2nd.nd$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.yau.2nd.nd$R2 < 0.9 | SR.yau.2nd.nd$p_value > 0.05
SR.yau.2nd.nd$`Sampling_Rate (m3/d)`[mask] <- NA
SR.yau.2nd.nd$R2[mask] <- NA
SR.yau.2nd.nd$p_value[mask] <- NA

# Export results
write.csv(SR.yau.2nd.nd,
          file = "Output/Data/csv/SamplingRates/Personal/SR.yau.2nd.nd.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.yau.2nd.nd$congener <- factor(SR.yau.2nd.nd$congener,
                              levels = unique(SR.yau.2nd.nd$congener))

ggplot(SR.yau.2nd.nd, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Ya'u SR vs logKoa regression 1 ------------------------------------------
# (1) Average both d and nd
sr.ave.yau <- as.data.frame(rowMeans(cbind(SR.yau.1st.nd$`Sampling_Rate (m3/d)`, 
                                              SR.yau.1st.nd$`Sampling_Rate (m3/d)`), 
                                        na.rm = TRUE))

sr.ave.yau$logkoa <- logKoa$logKoa
colnames(sr.ave.yau) <- c('ave_sr', 'logKoa')
# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.1 <- lm(log(sr.ave.yau$ave_sr) ~ sr.ave.yau$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.1)[1])  # exponentiate the intercept
b <- coef(model.yau.1)[2]       # coefficient for logKoa
r2 <- summary(model.yau.1)$r.squared

# plot
p.sr.yau.koa.1 <- ggplot(sr.ave.yau, aes(x = logKoa, y = ave_sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = 7.3, y = 7,
           label = paste("Vol. 3 (nd 1st & 2nd weeks)"),size = 5) +
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

p.sr.yau.koa.1

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.nd.1.png",
       plot = p.sr.yau.koa.1, width = 6, height = 6, dpi = 500)

# (2) Individual values
# Create a long dataframe combining SR.yau.1st. nd and SR.yau.2nd.nd
sr.long.yau <- data.frame(
  sr = c(SR.yau.1st.nd$`Sampling_Rate (m3/d)`, SR.yau.2nd.nd$`Sampling_Rate (m3/d)`),
  logKoa = rep(logKoa$logKoa, 2)  # Repeat logKoa values for both d and nd
)

# Remove any NA values
sr.long.yau <- na.omit(sr.long.yau)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.2 <- lm(log(sr.long.yau$sr) ~ sr.long.yau$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.2)[1])  # exponentiate the intercept
b <- coef(model.yau.2)[2]       # coefficient for logKoa
r2 <- summary(model.yau.2)$r.squared

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

# Plot
p.sr.yau.koa.2 <- ggplot(sr.long.yau, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.long.yau$logKoa) + 1.2, y = max(sr.long.yau$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.long.yau$logKoa) + 0.35, y = max(sr.long.yau$sr) * 1.13,
           label = paste("R² =", round(r2, 2)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.yau.koa.2

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.nd.2.png",
       plot = p.sr.yau.koa.2, width = 6, height = 6, dpi = 500)

# (3) 1st week only
sr.long.yau.1st <- data.frame(
  sr = c(SR.yau.1st.nd$`Sampling_Rate (m3/d)`), logKoa = logKoa$logKoa)

# Remove any NA values
sr.long.yau.1st <- na.omit(sr.long.yau.1st)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.3 <- lm(log(sr.long.yau.1st$sr) ~ sr.long.yau.1st$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.3)[1])  # exponentiate the intercept
b <- coef(model.yau.3)[2]       # coefficient for logKoa
r2 <- summary(model.yau.3)$r.squared

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

# Plot
p.sr.yau.koa.3 <- ggplot(sr.long.yau.1st, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  xlim(7, 10) +
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

p.sr.yau.koa.3

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.nd.1st.week.png",
       plot = p.sr.yau.koa.3, width = 6, height = 6, dpi = 500)

# (4) 2nd week only
sr.long.yau.2nd <- data.frame(
  sr = c(SR.yau.2nd.nd$`Sampling_Rate (m3/d)`), logKoa = logKoa$logKoa)

# Remove any NA values
sr.long.yau.2nd <- na.omit(sr.long.yau.2nd)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.yau.4 <- lm(log(sr.long.yau.2nd$sr) ~ sr.long.yau.2nd$logKoa)

# Get the coefficients
a <- exp(coef(model.yau.4)[1])  # exponentiate the intercept
b <- coef(model.yau.4)[2]       # coefficient for logKoa
r2 <- summary(model.yau.4)$r.squared

# Print equation
cat("Exponential Equation: sr = ", round(a, 3), " * exp(", round(b, 2), " * logKoa)\n")
cat("R² = ", round(r2, 2), "\n")

# Plot
p.sr.yau.koa.4 <- ggplot(sr.long.yau.2nd, aes(x = logKoa, y = sr)) +
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

p.sr.yau.koa.4

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/Yau_logKoa.nd.2nd.week.png",
       plot = p.sr.yau.koa.4, width = 6, height = 6, dpi = 500)

# Calculate personal sampling rate Ya'u 2nd -------------------------------
# 3 WBs were not wiped and 3 WBs were wiped
# Dominant hand used here
# WBs were used to calculate PCB concentration
# Concentrations were calculated for each sampling day 
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.yau.2 <- data.yau2[1:3, 5:177]
  # Calculate air concentration in ng/m3
  time <- data.yau2[1:3, 1]
  conc <- t(sweep(data.yau.2, 1, 0.5 * time, "/"))
  
  # Calculate Veff for WB nw
  mass.WD.nw <- data.yau2[4:6, 5:177]
  Veff.yau.nw.d <- mass.WD.nw/t(conc)
  # Add metadata to Veff and change format
  Veff.yau.nw.d <- cbind(data.yau2[4:6, 4], data.yau2[4:6, 1], Veff.yau.nw.d)
  # Transform to data.frame
  Veff.yau.nw.d <- as.data.frame(Veff.yau.nw.d)
  # Add names to first 2 columns
  colnames(Veff.yau.nw.d)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.yau.nw.d[, 2:175] <- apply(Veff.yau.nw.d[, 2:175], 2, as.numeric)
  
  # Calculate Veff for WB w
  mass.WD.w <- data.yau2[7:9, 5:177]
  Veff.yau.w <- mass.WD.w/t(conc)
  # Add metadata to Veff and change format
  Veff.yau.w.d <- cbind(data.yau2[7:9, 4], data.yau2[7:9, 1], Veff.yau.w)
  # Transform to data.frame
  Veff.yau.w.d <- as.data.frame(Veff.yau.w.d)
  # Add names to first 2 columns
  colnames(Veff.yau.w.d)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.yau.w.d[, 2:175] <- apply(Veff.yau.w.d[, 2:175], 2, as.numeric)
  # Select time
  Veff.yau.nw.d.t <- Veff.yau.nw.d[, 2]
  # Select time
  Veff.yau.w.d.t <- Veff.yau.w.d[, 2]
}

# Calculate sampling rate (SR) for nw & w (m3/d)
# (1) Remove metadata from Veff.yau.nw
Veff.yau.nw.d.2 <- Veff.yau.nw.d[, 3:175]
# Create matrix for sampling rate (SR)
SR.yau.nw.d <- matrix(nrow = length(Veff.yau.nw.d.2[1,]), ncol = 3)

for(i in 1:length(SR.yau.nw.d[, 1])) {
  if (sum(!is.na(Veff.yau.nw.d.2[, i]) & !is.infinite(Veff.yau.nw.d.2[, i])) == 3) {
    fit <- lm(Veff.yau.nw.d.2[, i] ~ 0 + Veff.yau.nw.d.t)
    SR.yau.nw.d[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.nw.d[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.nw.d[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.nw.d[i, 1] <- 0
    SR.yau.nw.d[i, 2] <- 0
    SR.yau.nw.d[i, 3] <- 0
  }
}

SR.yau.nw.d <- data.frame(SR.yau.nw.d, group = "ParticipantY.nw.d")
colnames(SR.yau.nw.d) <-c("Sampling_Rate (m3/d)", "R2", "p_value", "group")
congener <- names(head(Veff.yau.nw.d.2)[0, ])
SR.yau.nw.d <- cbind(congener, SR.yau.nw.d)

# Convert R2 and p-value to numeric
SR.yau.nw.d$`Sampling_Rate (m3/d)` <- as.numeric(SR.yau.nw.d$`Sampling_Rate (m3/d)`)
SR.yau.nw.d$R2 <- as.numeric(SR.yau.nw.d$R2)
SR.yau.nw.d$p_value <- as.numeric(SR.yau.nw.d$p_value)

# Update R2 and p-value to NA based on conditions
mask <- SR.yau.nw.d$R2 < 0.9 | SR.yau.nw.d$p_value > 0.05
SR.yau.nw.d$`Sampling_Rate (m3/d)`[mask] <- NA
SR.yau.nw.d$R2[mask] <- NA
SR.yau.nw.d$p_value[mask] <- NA

# Export results
write.csv(SR.yau.nw.d,
          file = "Output/Data/csv/SamplingRates/Personal/SR.yau.nw.d.csv",
          row.names = FALSE)

# Plot
# Organize PCB names
SR.yau.nw.d$congener <- factor(SR.yau.nw.d$congener,
                              levels = unique(SR.yau.nw.d$congener))

ggplot(SR.yau.nw.d, aes(x = congener, y = `Sampling_Rate (m3/d)`, color = group)) +
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

