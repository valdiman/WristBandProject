## Script to calcualte individual PCB sampling rates
# for silicone wristbands. Sampling rates were calculated with
# (i) static (2 times) and (ii) rotating set up.
# Airborne concentration was measured using low volume samplers (PUF)

# Install packages
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("dplyr")

# Load libraries
{
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
}

# Sampling rates calculations under static conditions ---------------------
# Read data
{
  PUF <- read.csv("Data/PUF.csv") # ng/m3
  WB <- read.csv("Data/WB.csv")
  logKoa <- read.csv("Data/logKoa.csv")
}

# Remove metadata
WB.1 <- subset(WB, select = -c(sample:time))

# Perform linear regression of WD accumulate mass vs time to check linearity
# Create matrix to storage data
WBMatrix <- matrix(nrow = length(WB.1), ncol = 3)

for(i in 1:length(WB.1)) {
  fit <- lm(WB.1[,i] ~ WB$time)
  WBMatrix[i,1] <- format(signif(summary(fit)$coef[2,"Estimate"], digits = 2))
  WBMatrix[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 4))
  WBMatrix[i,3] <- format(signif(summary(fit)$coef[2,"Pr(>|t|)"], digits = 3))
}

colnames(WBMatrix) <-c("slope", "R2", "p-value")
congener <- names(head(WB.1)[0,])
WBMatrix <- cbind(congener, WBMatrix)

# Export
write.csv(WBMatrix, file = "Output/Data/csv/SamplingRates/SR/WDLinearity.csv")

# Calculate sampling rate
# Remove sample name of PUFs
PUF.1 <- subset(PUF, select = -c(sample))
# Average individual PCBs from PUF
PUF.mean <- t(colMeans(PUF.1))
# Convert time into days
time <- WB$time/24

# (dMWD/Cair) = Rsdt
# Force intercept 0
# Conditions: Mean of PUF >0, but also all the measurements >0,
# Measurements from WB >3
# Create matrix for sampling rate (SR)
SR.st.1 <- matrix(nrow = length(WB.1), ncol = 3)

for(i in 1:length(WB.1)) {
  if ((PUF.mean[i] > 0) && (sum(WB.1[i] > 0,
                                na.rm = TRUE) > 3) && (sum(PUF.1[i] > 0,
                                                           na.rm = TRUE) == 4)) {
    fit <- lm(WB.1[,i]/PUF.mean[i] ~ 0 + time)
    SR.st.1[i,1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.st.1[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.st.1[i,3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.st.1[i,1] <- NA
    SR.st.1[i,2] <- NA
    SR.st.1[i,3] <- NA
  }
}

colnames(SR.st.1) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.1)[0,])
SR.st.1 <- as.data.frame(cbind(congener, SR.st.1))

# Convert R2 and p-value to numeric
SR.st.1$`Sampling Rate (m3/d)` <- as.numeric(SR.st.1$`Sampling Rate (m3/d)`)
SR.st.1$R2 <- as.numeric(SR.st.1$R2)
SR.st.1$`p-value` <- as.numeric(SR.st.1$`p-value`)

# Update R2 and p-value to NA based on conditions
mask <- SR.st.1$R2 < 0.9 | SR.st.1$`p-value` > 0.05
SR.st.1$`Sampling Rate (m3/d)`[mask] <- NA
SR.st.1$R2[mask] <- NA
SR.st.1$`p-value`[mask] <- NA
Awb <- 0.0054773 # [m2]
SR.st.1$ko <- SR.st.1$`Sampling Rate (m3/d)` /  Awb # [m/d]

# Values
SR.n <- length(na.omit(SR.st.1$`Sampling Rate (m3/d)`))
SR.ave <- mean(SR.st.1$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.sd <- sd(SR.st.1$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.ko.ave <- mean(SR.st.1$ko, na.rm = TRUE) 
SR.ko.sd <- sd(SR.st.1$ko, na.rm = TRUE)

# Add average of ko to the Na values.
SR.st.1$ko[is.na(SR.st.1$ko)] <- SR.ko.ave

# Export
write.csv(SR.st.1, file = "Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV1.csv")

# Plots
# Change number congener in [] 

fit1 <- lm(WB.1$PCB18.30/PUF.mean[16] ~ 0 + time)
ggplot(WB, aes(x = time, y = WB.1$PCB18.30/PUF.mean[16])) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  xlab(expression(bold("Deployment time (hr)"))) +
  ylab(expression(bold("Effective Volume PCBs 18+30 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 10)) +
  geom_point(color = "black", shape = 21, fill = "grey", size = 1.5) +
  stat_smooth(method = "lm", col = "black", se = FALSE,
              formula = y ~ 0 + x, fullrange = T) +
  xlim(0, 125) +
  ylim(0, 4) +
  annotate("text", x = 72, y = 0.2,
           label = paste(" Slope (m3/d)=",
                         signif(summary(fit1)$coef[1,"Estimate"], 2),
                         ", R2 = ", signif(summary(fit1)$adj.r.squared, 3)),
                         size = 2.5, fontface = 2)

# Sampling rates calculations under static and rotating conditions -------------------
# Read data
{
  PUF.v2 <- read.csv("Data/PUF2.csv") # ng/m3
  WB.st <- read.csv("Data/WB2.csv")
  WB.rot <- read.csv("Data/WBROT.csv")
}

# Remove metadata
WB.st.2 <- subset(WB.st, select = -c(sample:time))
WB.rot.2 <- subset(WB.rot, select = -c(sample:time))
# Convert time into days
time.2 <- WB.st$time/24

# Calculate sampling rate for stationary 2nd time
# Remove sample name of PUFs
PUF.2 <- subset(PUF.v2, select = -c(Concentration..ng.m3.))
# Average individual PCBs from PUF
# If value 0, the average is the number left.
PUF.mean.2 <- t(apply(PUF.2, 2, function(x) mean(x[x != 0])))
colnames(PUF.mean.2) <- colnames(PUF.2)  # Ensure the column names match the original data

# (dMWD/Cair) = Rsdt
# Force intercept 0
# Conditions: Mean of PUF >0, but also all the measurements >0,
# Measurements from WB >3
# Create matrix for sampling rate (SR)
SR.st.2 <- matrix(nrow = length(WB.st.2), ncol = 3)

for(i in 1:length(WB.st.2)) {
  if (!is.na(PUF.mean.2[i]) && PUF.mean.2[i] > 0) { # for NaN values
    fit <- lm(WB.st.2[,i]/PUF.mean.2[i] ~ 0 + time.2)
    SR.st.2[i,1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.st.2[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.st.2[i,3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.st.2[i,1] <- NA
    SR.st.2[i,2] <- NA
    SR.st.2[i,3] <- NA
  }
}

colnames(SR.st.2) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.st.2)[0,])
SR.st.2 <- as.data.frame(cbind(congener, SR.st.2))

# Convert R2 and p-value to numeric
SR.st.2$`Sampling Rate (m3/d)` <- as.numeric(SR.st.2$`Sampling Rate (m3/d)`)
SR.st.2$R2 <- as.numeric(SR.st.2$R2)
SR.st.2$`p-value` <- as.numeric(SR.st.2$`p-value`)

# Update R2 and p-value to NA based on conditions
mask <- SR.st.2$R2 < 0.9 | SR.st.2$`p-value` > 0.05
SR.st.2$`Sampling Rate (m3/d)`[mask] <- NA
SR.st.2$R2[mask] <- NA
SR.st.2$`p-value`[mask] <- NA
SR.st.2$ko <- SR.st.2$`Sampling Rate (m3/d)` / Awb # [m/d].

# Values
SR.st.n <- length(na.omit(SR.st.2$`Sampling Rate (m3/d)`))
SR.st.ave <- mean(SR.st.2$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.st.sd <- sd(SR.st.2$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.ko.ave <- mean(SR.st.2$ko, na.rm = TRUE) 
SR.ko.sd <- sd(SR.st.2$ko, na.rm = TRUE)

# Add average of ko to the Na values.
SR.st.2$ko[is.na(SR.st.2$ko)] <- SR.ko.ave

#export
write.csv(SR.st.2, file = "Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV2.csv")

# Plots
# Change number congener in [] 

fit1 <- lm(WB.st.2[, 153]/PUF.mean.2[153] ~ 0 + time.2)
ggplot(WB.st.2, aes(x = time.2*24, y = WB.st.2[, 153]/PUF.mean.2[153])) +
  theme_bw() +
  theme(aspect.ratio = 10/15) +
  xlab(expression(bold("Deployment time (hr)"))) +
  ylab(expression(bold("Effective Volume PCB 187 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 10)) +
  geom_point(color = "black", shape = 21, fill = "grey", size = 1.5) +
  stat_smooth(method = "lm", col = "black", se = FALSE,
              formula = y ~ 0 + x, fullrange = T) +
  xlim(0, 125) +
  ylim(0, 3) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 70, y = 3,
           label = paste(" Slope (m3/d)=",
                         signif(summary(fit1)$coef[1,"Estimate"], 2),
                         ", R2 = ", signif(summary(fit1)$adj.r.squared, 3)),
           size = 2.5, fontface = 2)

# Calculate sampling rate for rotating
# (dMWD/Cair) = Rsdt
# Force intercept 0
# Conditions: Mean of PUF >0, but also all the measurements > 0,
# Measurements from WB >3
# Create matrix for sampling rate (SR)
SR.rot <- matrix(nrow = length(WB.rot.2), ncol = 3)

for(i in 1:length(WB.rot.2)) {
  if (!is.na(PUF.mean.2[i]) && PUF.mean.2[i] > 0) { # for NaN values
    fit <- lm(WB.rot.2[,i]/PUF.mean.2[i] ~ 0 + time.2)
    SR.rot[i,1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.rot[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.rot[i,3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.rot[i,1] <- NA
    SR.rot[i,2] <- NA
    SR.rot[i,3] <- NA
  }
}

colnames(SR.rot) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.rot.2)[0,])
SR.rot <- as.data.frame(cbind(congener, SR.rot))

# Convert R2 and p-value to numeric
SR.rot$`Sampling Rate (m3/d)` <- as.numeric(SR.rot$`Sampling Rate (m3/d)`)
SR.rot$R2 <- as.numeric(SR.rot$R2)
SR.rot$`p-value` <- as.numeric(SR.rot$`p-value`)

# Update R2 and p-value to NA based on conditions
mask <- SR.rot$R2 < 0.9 | SR.rot$`p-value` > 0.05
SR.rot$`Sampling Rate (m3/d)`[mask] <- NA
SR.rot$R2[mask] <- NA
SR.rot$`p-value`[mask] <- NA
SR.rot$ko <- SR.rot$`Sampling Rate (m3/d)` / Awb # [m/d]

# Values
SR.rot.n <- length(na.omit(SR.rot$`Sampling Rate (m3/d)`))
SR.rot.ave <- mean(SR.rot$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.rot.sd <- sd(SR.rot$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.ko.ave <- mean(SR.rot$ko, na.rm = TRUE) 
SR.ko.sd <- sd(SR.rot$ko, na.rm = TRUE)

# Add averege of ko to the Na values.
SR.rot$ko[is.na(SR.rot$ko)] <- SR.ko.ave

# Export
write.csv(SR.rot, file = "Output/Data/csv/SamplingRates/SR/WDSamplingRateRotV1.csv")

# Plots
# Change number congener in [] 

fit1 <- lm(WB.rot.2$PCB18.30/PUF.mean.2[17] ~ 0 + time.2)
ggplot(WB.rot.2, aes(x = time.2*24, y = `PCB18.30`/PUF.mean.2[17])) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  xlab(expression(bold("Deployment time (hr)"))) +
  ylab(expression(bold("Effective Volume PCBs 18 + 30 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 10)) +
  geom_point(color = "black", shape = 21, fill = "grey", size = 1.5) +
  stat_smooth(method = "lm", col = "black", se = FALSE,
              formula = y ~ 0 + x, fullrange = T) +
  xlim(0, 125) +
  ylim(0, 6) +
  theme(plot.title = element_text(size = 8, face = "bold")) +
  annotate("text", x = 50, y = 5.5,
           label = paste(" Slope (m3/d)=",
                         signif(summary(fit1)$coef[1,"Estimate"], 2),
                         ", R2 = ", signif(summary(fit1)$adj.r.squared, 3)),
           size = 2.5, fontface = 2)

# Plot both 'stat' and 'dyn' sampling rates -------------------------------
# Combine plot
# PCB 18+30
WB$PCB18.30_comb <- WB.1$PCB18.30 / PUF.mean[16]
WB.rot.2$PCB18.30_comb <- WB.rot.2$PCB18.30 / PUF.mean.2[17]

# Combine data, padding shorter dataset with NA
max_len <- max(nrow(WB), nrow(WB.rot.2))

WB_padded <- data.frame(
  time = c(WB$time, rep(NA, max_len - nrow(WB))),
  PCB18.30_comb = c(WB$PCB18.30_comb, rep(NA, max_len - nrow(WB))),
  group = rep("Static", max_len)
)

WB_rot_padded <- data.frame(
  time = c(WB.rot$time, rep(NA, max_len - nrow(WB.rot.2))),
  PCB18.30_comb = c(WB.rot.2$PCB18.30_comb, rep(NA, max_len - nrow(WB.rot.2))),
  group = rep("Dynamic", max_len)
)

combined_data <- rbind(WB_padded, WB_rot_padded)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(PCB18.30_comb ~ 0 + time))[1] * 24) # times 24 to get [m3/d]

# Plotting
plot.18.30 <- ggplot(combined_data, aes(x = time, y = PCB18.30_comb,
                                        color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  annotate("text", x = 1, y = 5.5,
           label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0,
           label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  theme_bw() +
  xlim(0, 125) +
  ylim(0, 6) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCBs 18+30 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        aspect.ratio = 1.5)

# See plot
plot.18.30

# Save plot
ggsave("Output/Plots/SamplingRates/SR/PCB18.30SamplingRates.png",
       plot = plot.18.30, width = 8, height = 10, dpi = 1300)

# PCB 52
# This plot is included in the paper (Fig. 2)
WB$PCB52_comb <- WB.1$PCB52 / PUF.mean[45]
WB.rot.2$PCB52_comb <- WB.rot.2$PCB52 / PUF.mean.2[45]

# Combine data, padding shorter dataset with NA
max_len <- max(nrow(WB), nrow(WB.rot.2))

WB_padded <- data.frame(
  time = c(WB$time, rep(NA, max_len - nrow(WB))),
  PCB52_comb = c(WB$PCB52_comb, rep(NA, max_len - nrow(WB))),
  group = rep("Static", max_len)
)

WB_rot_padded <- data.frame(
  time = c(WB.rot$time, rep(NA, max_len - nrow(WB.rot.2))),
  PCB52_comb = c(WB.rot.2$PCB52_comb, rep(NA, max_len - nrow(WB.rot.2))),
  group = rep("Dynamic", max_len)
)

combined_data <- rbind(WB_padded, WB_rot_padded)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(PCB52_comb ~ 0 + time))[1] * 24) # times 24 to get [m3/d]

# Define colors for groups (Dynamic and Static)
group_colors <- c("Dynamic" = "#F8766D",
                  "Static" = "#00BFC4")

# Plot with defined colors
plot.52 <- ggplot(combined_data, aes(x = time, y = PCB52_comb,
                                     color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  theme_bw() +
  xlim(0, 125) +
  ylim(0, 6) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 52 (m"^3*")"))) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.y = element_text(face = "bold", size = 28),
        axis.title.y = element_text(face = "bold", size = 28),
        axis.text.x = element_text(face = "bold", size = 28),
        axis.title.x = element_text(face = "bold", size = 28),
        legend.position = "none",  # Remove the existing legend
        aspect.ratio = 1.5)

# Adding color symbols before text with increased spacing
plot.52 <- plot.52 + 
  # Add color dot for "Dynamic"
  geom_point(aes(x = 1, y = 6), size = 6, shape = 21,
             fill = group_colors["Dynamic"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 5, y = 6,
           label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  # Add color dot for "Static"
  geom_point(aes(x = 1, y = 5.5), size = 6, shape = 21,
             fill = group_colors["Static"], color = "black",
             show.legend = FALSE) +
  annotate("text", x = 5, y = 5.5,
           label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  # Add label for "Study 1"
  annotate("text", x = 90, y = 0.1,
           label = "Study 1", hjust = 0, size = 10, color = "black") +
  # Add label for panel (a)
  annotate("text", x = Inf, y = Inf,
           label = "(a)", hjust = 1.1, vjust = 1.5, 
           size = 10, color = "black")

# See plot
plot.52

# Save plot
ggsave("Output/Plots/SamplingRates/SR/PCB52SamplingRates5.png",
       plot = plot.52, width = 8, height = 10, dpi = 1300)

# PCB 118
WB$PCB118_comb <- WB.1$PCB118 / PUF.mean[96]
WB.rot.2$PCB118_comb <- WB.rot.2$PCB118 / PUF.mean.2[96]

# Combine data, padding shorter dataset with NA
max_len <- max(nrow(WB), nrow(WB.rot.2))

WB_padded <- data.frame(
  time = c(WB$time, rep(NA, max_len - nrow(WB))),
  PCB118_comb = c(WB$PCB118_comb, rep(NA, max_len - nrow(WB))),
  group = rep("Static", max_len)
)

WB_rot_padded <- data.frame(
  time = c(WB.rot$time, rep(NA, max_len - nrow(WB.rot.2))),
  PCB118_comb = c(WB.rot.2$PCB118_comb, rep(NA, max_len - nrow(WB.rot.2))),
  group = rep("Dynamic", max_len)
)

combined_data <- rbind(WB_padded, WB_rot_padded)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(PCB118_comb ~ 0 + time))[1] * 24) # times 24 to get [m3/d]

# Plotting
plot.118 <- ggplot(combined_data, aes(x = time, y = PCB118_comb,
                                      color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  annotate("text", x = 1, y = 5.5,
           label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0,
           label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  theme_bw() +
  xlim(0, 125) +
  ylim(0, 6) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 118 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        aspect.ratio = 1.5)

# See plot
plot.118

# Save plot
ggsave("Output/Plots/SamplingRates/SR/PCB118SamplingRates.png",
       plot = plot.118, width = 8, height = 10, dpi = 1300)

# PCB 187
WB$PCB187_comb <- WB.1$PCB187 / PUF.mean[152]
WB.rot.2$PCB187_comb <- WB.rot.2$PCB187 / PUF.mean.2[153]

# Combine data, padding shorter dataset with NA
max_len <- max(nrow(WB), nrow(WB.rot.2))

WB_padded <- data.frame(
  time = c(WB$time, rep(NA, max_len - nrow(WB))),
  PCB187_comb = c(WB$PCB187_comb, rep(NA, max_len - nrow(WB))),
  group = rep("Static", max_len)
)

WB_rot_padded <- data.frame(
  time = c(WB.rot$time, rep(NA, max_len - nrow(WB.rot.2))),
  PCB187_comb = c(WB.rot.2$PCB187_comb, rep(NA, max_len - nrow(WB.rot.2))),
  group = rep("Dynamic", max_len)
)

combined_data <- rbind(WB_padded, WB_rot_padded)

slopes <- combined_data %>%
  group_by(group) %>%
  summarize(slope = coef(lm(PCB187_comb ~ 0 + time))[1] * 24) # times 24 to get [m3/d]

# Plotting
plot.187 <- ggplot(combined_data, aes(x = time, y = PCB187_comb,
                                      color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  annotate("text", x = 1, y = 5.5,
           label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0,
           label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  theme_bw() +
  xlim(0, 125) +
  ylim(0, 6) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 187 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        aspect.ratio = 1.5)

# See plot
plot.187

# Save plot
ggsave("Output/Plots/SamplingRates/SR/PCB187SamplingRates.png",
       plot = plot.187, width = 8, height = 10, dpi = 1300)

# SR vs logKoa regressions ------------------------------------------------
# (1) Static 1
# Remove PCB14 & PCBs128+166
logKoa.0 <- logKoa[!logKoa$congener %in% c("PCB14", "PCB128.166"), ]
sr.stat.1 <- data.frame(
  sr = SR.st.1$`Sampling Rate (m3/d)`,
  logKoa = logKoa.0$logKoa)

# Remove any NA values
sr.stat.1 <- na.omit(sr.stat.1)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.stat.1 <- lm(log(sr.stat.1$sr) ~ sr.stat.1$logKoa)

# Get the coefficients
a <- exp(coef(model.stat.1)[1])  # exponentiate the intercept
b <- coef(model.stat.1)[2]       # coefficient for logKoa
r2 <- summary(model.stat.1)$r.squared

# Plot
p.sr.stat.1.koa <- ggplot(sr.stat.1, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.stat.1$logKoa) + 1.05, y = max(sr.stat.1$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.stat.1$logKoa) + 0.35, y = max(sr.stat.1$sr) * 1.13,
           label = paste("R² =", round(r2, 3)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.stat.1.koa

# Save plot in folder
ggsave("Output/Plots/SamplingRates/SR/Stat1_logKoa.png", plot = p.sr.stat.1.koa,
       width = 6, height = 6, dpi = 500)

# (2) Static 2
logKoa.1 <- logKoa[logKoa$congener != "PCB14", ]
sr.stat.2 <- data.frame(
  sr = SR.st.2$`Sampling Rate (m3/d)`,
  logKoa = logKoa.1$logKoa)

# Remove any NA values
sr.stat.2 <- na.omit(sr.stat.2)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.stat.2 <- lm(log(sr.stat.2$sr) ~ sr.stat.2$logKoa)

# Get the coefficients
a <- exp(coef(model.stat.2)[1])  # exponentiate the intercept
b <- coef(model.stat.2)[2]       # coefficient for logKoa
r2 <- summary(model.stat.2)$r.squared

# Plot
p.sr.stat.2.koa <- ggplot(sr.stat.2, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.stat.2$logKoa) + 1.05, y = max(sr.stat.2$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 2), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.stat.2$logKoa) + 0.35, y = max(sr.stat.2$sr) * 1.13,
           label = paste("R² =", round(r2, 3)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.stat.2.koa

# Save plot in folder
ggsave("Output/Plots/SamplingRates/SR/Stat2_logKoa.png", plot = p.sr.stat.2.koa,
       width = 6, height = 6, dpi = 500)

# (3) Rot
sr.rot <- data.frame(
  sr = SR.rot$`Sampling Rate (m3/d)`,
  logKoa = logKoa.1$logKoa)

# Remove any NA values
sr.rot <- na.omit(sr.rot)

# Fit exponential regression model: sr = a * exp(b * logKoa)
model.rot <- lm(log(sr.rot$sr) ~ sr.rot$logKoa)

# Get the coefficients
a <- exp(coef(model.rot)[1])  # exponentiate the intercept
b <- coef(model.rot)[2]       # coefficient for logKoa
r2 <- summary(model.rot)$r.squared

# Plot
p.sr.rot.koa <- ggplot(sr.rot, aes(x = logKoa, y = sr)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = FALSE, color = "blue") +
  annotate("text", x = min(sr.rot$logKoa) + 1.05, y = max(sr.rot$sr) * 1.2,
           label = paste("sr =", round(a, 3), "* exp(", round(b, 3), "* log Koa)"),
           size = 5) +
  annotate("text", x = min(sr.rot$logKoa) + 0.35, y = max(sr.rot$sr) * 1.13,
           label = paste("R² =", round(r2, 5)), size = 5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab(expression(bold("log Koa"))) +
  ylab(expression(bold("Ave Sampling Rate (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

p.sr.rot.koa

# Save plot in folder
ggsave("Output/Plots/SamplingRates/SR/Rot_logKoa.png", plot = p.sr.rot.koa,
       width = 6, height = 6, dpi = 500)
