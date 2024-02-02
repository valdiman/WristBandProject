## Code to calculate individual PCB "personal" sampling rates
# for silicone wristbands.

# Install packages
install.packages("readxl") #say no!
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

# Calculate personal sampling rate Amanda ---------------------------------
# WBs were used to calculate PCB concentration
# triplicates for 4.3 days were deployed
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.amanda.1 <- data.amanda[1:3,]
  # Average 3 WBs
  data.amanda.2 <- colMeans(data.amanda.1[, 3:175])
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
  Veff.amanda.r <- Veff.amanda[1:5, 3:175]
  # Select time
  Veff.amanda.r.t <- Veff.amanda[1:5, 2]
  # Select left, remove metadata
  Veff.amanda.l <- Veff.amanda[6:10, 3:175]
  # Select time
  Veff.amanda.l.t <- Veff.amanda[6:10, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.amanda.r <- matrix(nrow = length(Veff.amanda.r[1,]), ncol = 3)

for(i in 1:length(SR.amanda.r[, 1])) {
  if (length(unique(Veff.amanda.r[, i])) >= 3) {
    fit <- lm(Veff.amanda.r[, i] ~ 0 + Veff.amanda.r.t)
    SR.amanda.r[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.amanda.r[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.amanda.r[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.amanda.r[i, 1] <- 0
    SR.amanda.r[i, 2] <- 0
    SR.amanda.r[i, 3] <- 0
  }
}

colnames(SR.amanda.r) <-c("Sampling_Rate", "R2", "p_value")
congener <- names(head(Veff.amanda.r)[0, ])
SR.amanda.r <- cbind(congener, SR.amanda.r)
SR.amanda.r <- data.frame(SR.amanda.r, group = "ParticipantA.r")

# Create matrix for sampling rate (SR)
SR.amanda.l <- matrix(nrow = length(Veff.amanda.l[1,]), ncol = 3)

for(i in 1:length(SR.amanda.l[, 1])) {
  if (length(unique(Veff.amanda.l[, i])) >= 3) {
    fit <- lm(Veff.amanda.l[, i] ~ 0 + Veff.amanda.l.t)
    SR.amanda.l[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.amanda.l[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.amanda.l[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.amanda.l[i, 1] <- 0
    SR.amanda.l[i, 2] <- 0
    SR.amanda.l[i, 3] <- 0
  }
}

colnames(SR.amanda.l) <-c("Sampling_Rate", "R2", "p_value")
SR.amanda.l <- cbind(congener, SR.amanda.l)
SR.amanda.l <- data.frame(SR.amanda.l, group = "ParticipantA.l")

# Plot
# Convert numerical columns to numeric
SR.amanda.r[, 2:4] <- apply(SR.amanda.r[, 2:4], 2, as.numeric)
# Organize PCB names
SR.amanda.r$congener <- factor(SR.amanda.r$congener,
                            levels = unique(SR.amanda.r$congener))

# Plot with legend
ggplot(SR.amanda.r[SR.amanda.r$Sampling_Rate > 0 & SR.amanda.r$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Convert numerical columns to numeric
SR.amanda.l[, 2:4] <- apply(SR.amanda.l[, 2:4], 2, as.numeric)
# Organize PCB names
SR.amanda.l$congener <- factor(SR.amanda.l$congener,
                               levels = unique(SR.amanda.l$congener))

ggplot(SR.amanda.l[SR.amanda.l$Sampling_Rate > 0 & SR.amanda.l$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Calculate personal sampling rate Kay ------------------------------------
# WBs were used to calculate PCB concentration
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
  Veff.kay.r <- Veff.kay[1:5, 3:175]
  # Select time
  Veff.kay.r.t <- Veff.kay[1:5, 2]
}

# Calculate sampling rate (SR) for right and left hands (m3/d)
# Create matrix for sampling rate (SR)
SR.kay.r <- matrix(nrow = length(Veff.kay.r[1,]), ncol = 3)

for(i in 1:length(SR.kay.r[, 1])) {
  if (length(unique(Veff.kay.r[, i])) >= 3) {
    fit <- lm(Veff.kay.r[, i] ~ 0 + Veff.kay.r.t)
    SR.kay.r[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.kay.r[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.kay.r[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.kay.r[i, 1] <- 0
    SR.kay.r[i, 2] <- 0
    SR.kay.r[i, 3] <- 0
  }
}

colnames(SR.kay.r) <-c("Sampling_Rate", "R2", "p_value")
SR.kay.r <- cbind(congener, SR.kay.r)
SR.kay.r <- data.frame(SR.kay.r, group = "ParticipantK.r")

# Plot
SR.kay.r[, 2:4] <- apply(SR.kay.r[, 2:4], 2, as.numeric)
# Organize PCB names
SR.kay.r$congener <- factor(SR.kay.r$congener,
                               levels = unique(SR.kay.r$congener))

ggplot(SR.kay.r[SR.kay.r$Sampling_Rate > 0 & SR.kay.r$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Calculate personal sampling rate Ya'u -----------------------------------
# WBs were used to calculate PCB concentration
# Concentrations were calculated for each sampling day 
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.yau.1 <- data.yau[1:6, 4:174]
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  time <- data.yau[1:6, 1]
  conc <- t(sweep(data.yau.1, 1, 0.5 * time, "/"))
  # get WB mass
  mass.WD <- data.yau[7:12, 4:174]
  Veff.yau <- mass.WD/t(conc)
  # Add metadata to Veff.amanda and change format
  Veff.yau <- cbind(data.yau[7:12, 3], data.yau[7:12, 1], Veff.yau)
  # Transform to data.frame
  Veff.yau <- as.data.frame(Veff.yau)
  # Add names to first 2 columns
  colnames(Veff.yau)[1:2] <- c("sample", "time.day")
  # Change characters to numbers format
  Veff.yau[, 2:173] <- apply(Veff.yau[, 2:173], 2, as.numeric)
  # Select 1st week, remove metadata
  Veff.yau.1st <- Veff.yau[1:3, 3:173]
  # Select time
  Veff.yau.1st.t <- Veff.yau[1:3, 2]
  # Select left, remove metadata
  Veff.yau.2nd <- Veff.yau[4:6, 3:173]
  # Select time
  Veff.yau.2nd.t <- Veff.yau[4:6, 2]
}

# Calculate sampling rate (SR) for 1st and 2nd weeks (m3/d)
# Create matrix for sampling rate (SR)
SR.yau.1st <- matrix(nrow = length(Veff.yau.1st[1,]), ncol = 3)

for(i in 1:length(SR.yau.1st[, 1])) {
  if (sum(!is.na(Veff.yau.1st[, i ]) & !is.infinite(Veff.yau.1st[, i])) == 3) {
    fit <- lm(Veff.yau.1st[, i] ~ 0 + Veff.yau.1st.t)
    SR.yau.1st[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.1st[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.1st[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.1st[i, 1] <- 0
    SR.yau.1st[i, 2] <- 0
    SR.yau.1st[i, 3] <- 0
  }
}

colnames(SR.yau.1st) <-c("Sampling_Rate", "R2", "p_value")
congener <- names(head(Veff.yau.1st)[0, ])
SR.yau.1st <- cbind(congener, SR.yau.1st)
SR.yau.1st <- data.frame(SR.yau.1st, group = "ParticipantY.1st")

# Plot
SR.yau.1st[, 2:4] <- apply(SR.yau.1st[, 2:4], 2, as.numeric)
# Organize PCB names
SR.yau.1st$congener <- factor(SR.yau.1st$congener,
                            levels = unique(SR.yau.1st$congener))

ggplot(SR.yau.1st[SR.yau.1st$Sampling_Rate > 0 & SR.yau.1st$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
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
SR.yau.2nd <- matrix(nrow = length(Veff.yau.2nd[1,]), ncol = 3)

for(i in 1:length(SR.yau.2nd[, 1])) {
  if (sum(!is.na(Veff.yau.2nd[, i ]) & !is.infinite(Veff.yau.2nd[, i])) == 3) {
    fit <- lm(Veff.yau.2nd[, i] ~ 0 + Veff.yau.2nd.t)
    SR.yau.2nd[i, 1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.yau.2nd[i, 2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.yau.2nd[i, 3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.yau.2nd[i, 1] <- 0
    SR.yau.2nd[i, 2] <- 0
    SR.yau.2nd[i, 3] <- 0
  }
}

colnames(SR.yau.2nd) <-c("Sampling_Rate", "R2", "p_value")
congener <- names(head(Veff.yau.2nd)[0, ])
SR.yau.2nd <- cbind(congener, SR.yau.2nd)
SR.yau.2nd <- data.frame(SR.yau.2nd, group = "ParticipantY.2nd")

# Plot
SR.yau.2nd[, 2:4] <- apply(SR.yau.2nd[, 2:4], 2, as.numeric)
# Organize PCB names
SR.yau.2nd$congener <- factor(SR.yau.2nd$congener,
                              levels = unique(SR.yau.2nd$congener))

ggplot(SR.yau.2nd[SR.yau.2nd$Sampling_Rate > 0 & SR.yau.2nd$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Combine sampling rates --------------------------------------------------
# Combine the three data frames
combined_SR <- rbind(SR.amanda.r, SR.amanda.l, SR.kay.r, SR.yau.1st,
                     SR.yau.2nd)

# Plot the combined data with different colors for each group
Plot.SRs <- ggplot(combined_SR[combined_SR$Sampling_Rate > 0 & combined_SR$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point(size = 1) +
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
ggsave("Output/Plots/SRsV02.png",
       plot = Plot.SRs, width = 15, height = 5, dpi = 500)

# Box plot
Plot.SRs.boxplot <- ggplot(combined_SR[combined_SR$Sampling_Rate > 0 & combined_SR$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate)) +
  geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# See plot
print(Plot.SRs.boxplot)
# Save plot
ggsave("Output/Plots/SRsBoxplotV02.png",
       plot = Plot.SRs.boxplot, width = 15, height = 5, dpi = 500)

# Potential regressions ---------------------------------------------------
# Read logKoa
logKoa <- data.frame(read_excel("Data/logKoa.xlsx", sheet = "logKoa",
                                  col_names = TRUE, col_types = NULL))

# Reshape the dataset into a long format
logKoa_long <- logKoa %>%
  gather(key = "variable", value = "logKoa", -congener)

# Combine the datasets based on the 'congener' column
combined_data <- merge(combined_SR, logKoa_long, by = "congener")

# Analisys all data
# Perform exponential regression with the specified condition
exp_reg <- nls(Sampling_Rate ~ a * exp(b * logKoa), 
               data = subset(combined_data, Sampling_Rate > 0 & p_value < 0.05),
               start = list(a = 1, b = 1))

# Extract coefficients
a_value <- coef(exp_reg)[["a"]]
b_value <- coef(exp_reg)[["b"]]

# Calculate R-squared value
RSS <- sum(residuals(exp_reg)^2)
TSS <- sum((combined_data$Sampling_Rate - mean(combined_data$Sampling_Rate))^2)
rsquared <- 1 - RSS/TSS

# Errors
# Calculate Residual Standard Error (RSE)
RSE <- sqrt(RSS / (length(residuals(exp_reg)) - 2))
# Calculate Mean Absolute Error (MAE)
MAE <- mean(abs(residuals(exp_reg)))
# Calculate Mean Squared Error (MSE)
MSE <- mean(residuals(exp_reg)^2)
# Calculate Root Mean Squared Error (RMSE)
RMSE <- sqrt(MSE)

# Create a data frame for text annotations
annotation_df <- data.frame(
  equation = sprintf("SR = %.3f * exp(%.3f * logKoa)", a_value, b_value),
  rsquared = sprintf("R-squared = %.3f", rsquared),
  x = max(combined_data$logKoa),  # Set x to the maximum logKoa value
  y = max(combined_data$Sampling_Rate),  # Set y to the maximum Sampling Rate value
  group = "exp regression"  # Specify a group for the text annotations
)

# Plot the data with the filtered exponential regression line
Plot.exp.regr.all <- ggplot(subset(combined_data, Sampling_Rate > 0 & p_value < 0.05),
       aes(x = logKoa, y = Sampling_Rate, color = group)) +
  geom_point() +
  geom_line(data = data.frame(logKoa = seq(min(combined_data$logKoa),
                                           max(combined_data$logKoa),
                                           length.out = 100)),
            aes(x = logKoa, y = a_value * exp(b_value * logKoa)), 
            color = "red", linetype = "solid", linewidth = 1) +  # Add exponential regression line
  geom_text(data = annotation_df, aes(x, y, label = equation, group = group),
            hjust = 2, vjust = 1, size = 4) +
  geom_text(data = annotation_df, aes(x, y, label = rsquared, group = group),
            hjust = 3.5, vjust = 3, size = 4) +
  theme_bw() +
  theme(aspect.ratio = 10/10) +
  xlab(expression(bold("log Koa (PCBi)"))) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 14)) +
  scale_color_discrete(name = "Participants")

# See plot
print(Plot.exp.regr.all)
# Save plot
ggsave("Output/Plots/SRExpRegresionAllV02.png",
       plot = Plot.exp.regr.all, width = 8, height = 8, dpi = 500)

# Select participant with high SR in the high PCBs
group1 <- "ParticipantA.r"
group2 <- "ParticipantA.l"

# Filter the data for the selected groups
selected_data <- combined_data[combined_data$group %in% c(group1, group2), ]

# Perform exponential regression with the specified condition
exp_reg <- nls(Sampling_Rate ~ a * exp(b * logKoa), 
               data = subset(selected_data, Sampling_Rate > 0 & p_value < 0.05),
               start = list(a = 1, b = 1))

# Extract coefficients
a_value <- coef(exp_reg)[["a"]]
b_value <- coef(exp_reg)[["b"]]

# Calculate R-squared value
RSS <- sum(residuals(exp_reg)^2)
TSS <- sum((selected_data$Sampling_Rate - mean(selected_data$Sampling_Rate))^2)
rsquared <- 1 - RSS/TSS

# Create a data frame for text annotations
annotation_df <- data.frame(
  equation = sprintf("SR = %.3f * exp(%.3f * logKoa)", a_value, b_value),
  rsquared = sprintf("R-squared = %.3f", rsquared),
  x = max(selected_data$logKoa),  # Set x to the maximum logKoa value
  y = max(selected_data$Sampling_Rate),  # Set y to the maximum Sampling Rate value
  group = "exp regression"  # Specify a group for the text annotations
)

# Plot the data with the filtered exponential regression line
Plot.exp.regr <- ggplot(subset(selected_data, Sampling_Rate > 0 & p_value < 0.05),
       aes(x = logKoa, y = Sampling_Rate, color = group)) +
  geom_point() +
  geom_line(data = data.frame(logKoa = seq(min(selected_data$logKoa),
                                           max(selected_data$logKoa),
                                           length.out = 100)),
            aes(x = logKoa, y = a_value * exp(b_value * logKoa)), 
            color = "red", linetype = "solid", linewidth = 1) +  # Add exponential regression line
  geom_text(data = annotation_df, aes(x, y, label = equation, group = group),
            hjust = 2, vjust = 5, size = 4) +
  geom_text(data = annotation_df, aes(x, y, label = rsquared, group = group),
            hjust = 3.5, vjust = 7, size = 4) +
  theme_bw() +
  theme(aspect.ratio = 10/10) +
  xlab(expression(bold("log Koa (PCBi)"))) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 14)) +
  scale_color_discrete(name = "Participants")

# See plot
print(Plot.exp.regr)
# Save plot
ggsave("Output/Plots/SRExpRegresionV01.png",
       plot = Plot.exp.regr, width = 8, height = 8, dpi = 500)

# Profiles ----------------------------------------------------------------
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
colnames(profile.amanda) <- c("congeners", "Stat.day5", "ParticipantA.r.day5",
                              "ParticipantA.l.day5")
profile.amanda$congeners <- factor(profile.amanda$congeners,
                                   levels = unique(profile.amanda$congeners))

# Reshape the data frame to long format
profile_long <- profile.amanda %>%
  pivot_longer(cols = c(Stat.day5, ParticipantA.r.day5, ParticipantA.l.day5),
               names_to = "Participant",
               values_to = "Value")

# Define color
palette <- brewer.pal(3, "Set1")

# Profile plot
Plot.prof.amanda <- ggplot(profile_long, aes(x = congeners, y = Value,
                                             fill = Participant)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  xlab("") +
  ylim(0, 0.12) +
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

# Save plot
ggsave("Output/Plots/ProfAmandaV01.png",
       plot = Plot.prof.amanda, width = 15, height = 5, dpi = 500)



