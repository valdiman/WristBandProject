## Code to calculate individual PCB "personal" sampling rates
# for silicone wristbands.
# nd non-dominant hand
# d dominant hand

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
}

# Read data from excel ----------------------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
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

# Gas sampling rate calculations ------------------------------------------
cte <- 10^-8
fg <- 1/(10^(logKoa$logKoa)*cte + 1)
sr.gas.nd <- SR.amanda.nd$`Sampling_Rate (m3/d)`* fg 
sr.gas.d <- SR.amanda.d$`Sampling_Rate (m3/d)`* fg
sr <- data.frame(
  congener = SR.amanda.nd$congener,
  sr.nd = SR.amanda.nd$`Sampling_Rate (m3/d)`,
  sr.d = SR.amanda.d$`Sampling_Rate (m3/d)`,
  sr.gas.d, sr.gas.nd
)

# Reshape the data into long format
sr_long <- sr %>%
  pivot_longer(cols = c(sr.nd, sr.gas.nd, sr.d, sr.gas.d),
               names_to = "Type", values_to = "Value")

sr_long$congener <- factor(sr_long$congener, levels = unique(sr_long$congener))

sr_long <- na.omit(sr_long)

# Filter data to include only sr.nd and sr.gas.nd
sr_filtered_nd <- sr_long %>% 
  filter(Type %in% c("sr.nd", "sr.gas.nd"))

# Plot both sr.nd and sr.gas.nd in the same graph
plot.sr.nd <- ggplot(sr_filtered_nd, aes(x = congener, y = Value, fill = Type)) +
  geom_col(position = "dodge") +  # Dodge to separate bars
  theme_bw() +
  theme(aspect.ratio = 1/4) +
  xlab(expression(bold(" "))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  scale_fill_manual(values = c("sr.nd" = "blue", "sr.gas.nd" = "red"))

plot.sr.nd

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/sr.nd.png", plot = plot.sr.nd,
       width = 10, height = 5, dpi = 500)

# Filter data to include only sr.nd and sr.gas.nd
sr_filtered_d <- sr_long %>% 
  filter(Type %in% c("sr.d", "sr.gas.d"))

# Plot both sr.nd and sr.gas.nd in the same graph
plot.sr.d <- ggplot(sr_filtered_d, aes(x = congener, y = Value, fill = Type)) +
  geom_col(position = "dodge") +  # Dodge to separate bars
  theme_bw() +
  theme(aspect.ratio = 1/4) +
  xlab(expression(bold(" "))) +
  ylab(expression(bold("Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 10)) +
  scale_fill_manual(values = c("sr.d" = "blue", "sr.gas.d" = "red"))

plot.sr.d

# Save plot in folder
ggsave("Output/Plots/SamplingRates/Personal/sr.d.png", plot = plot.sr.d,
       width = 10, height = 5, dpi = 500)
