## Code to individual PCB calculate sampling rates
# for silicone wristbands. Sampling rates were calculated with
# (i) static (2 times) and (ii) rotating set up.

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("dplyr")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
}

# Sampling rates calculations under static conditions ---------------------
# Read data from excel
PUF <- data.frame(read_excel("Data/DataStatic.xlsx", sheet = "PUFv2",
                               col_names = TRUE, col_types = NULL))

WB <- data.frame(read_excel("Data/DataStatic.xlsx", sheet = "WDv2",
                                col_names = TRUE, col_types = NULL))

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
write.csv(WBMatrix, file = "Output/Data/csv/WDLinearity.csv")

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
SR <- matrix(nrow = length(WB.1), ncol = 3)

for(i in 1:length(WB.1)) {
  if ((PUF.mean[i] > 0) && (sum(WB.1[i] > 0,
                                na.rm = TRUE) > 3) && (sum(PUF.1[i] > 0,
                                                           na.rm = TRUE) == 4)) {
    fit <- lm(WB.1[,i]/PUF.mean[i] ~ 0 + time)
    SR[i,1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR[i,3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR[i,1] <- NA
    SR[i,2] <- NA
    SR[i,3] <- NA
  }
}

colnames(SR) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.1)[0,])
SR <- as.data.frame(cbind(congener, SR))

# Convert R2 and p-value to numeric
SR$`Sampling Rate (m3/d)` <- as.numeric(SR$`Sampling Rate (m3/d)`)
SR$R2 <- as.numeric(SR$R2)
SR$`p-value` <- as.numeric(SR$`p-value`)

# Update R2 and p-value to NA based on conditions
mask <- SR$R2 < 0.9 | SR$`p-value` > 0.05
SR$`Sampling Rate (m3/d)`[mask] <- NA
SR$R2[mask] <- NA
SR$`p-value`[mask] <- NA

# Values
SR.n <- length(na.omit(SR$`Sampling Rate (m3/d)`))
SR.ave <- mean(SR$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.sd <- sd(SR$`Sampling Rate (m3/d)`, na.rm = TRUE)

# Export
write.csv(SR, file = "Output/Data/csv/WDSamplingRateV1.csv")

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
# Read data from excel
PUF.v2 <- data.frame(read_excel("Data/DataV05.xlsx", sheet = "PUFv3",
                             col_names = TRUE, col_types = NULL))

WB.st <- data.frame(read_excel("Data/DataV05.xlsx", sheet = "WDST",
                            col_names = TRUE, col_types = NULL))

WB.rot <- data.frame(read_excel("Data/DataV05.xlsx", sheet = "WDROT",
                               col_names = TRUE, col_types = NULL))

# Remove metadata
WB.st.2 <- subset(WB.st, select = -c(sample:time))
WB.rot.2 <- subset(WB.rot, select = -c(sample:time))
# Convert time into days
time.2 <- WB.st$time/24

# Calculate sampling rate for stationary 2nd time
# Remove sample name of PUFs
PUF.2 <- subset(PUF.v2, select = -c(Concentration..ng.m3.))
# Average individual PCBs from PUF
PUF.mean.2 <- t(colMeans(PUF.2))

# (dMWD/Cair) = Rsdt
# Force intercept 0
# Conditions: Mean of PUF >0, but also all the measurements >0,
# Measurements from WB >3
# Create matrix for sampling rate (SR)
SR.st <- matrix(nrow = length(WB.st.2), ncol = 3)

for(i in 1:length(WB.st.2)) {
  if (PUF.mean.2[i] > 0) {
    fit <- lm(WB.st.2[,i]/PUF.mean.2[i] ~ 0 + time.2)
    SR.st[i,1] <- format(signif(summary(fit)$coef[1,"Estimate"], digits = 3))
    SR.st[i,2] <- format(signif(summary(fit)$adj.r.squared, digits = 3))
    SR.st[i,3] <- format(signif(summary(fit)$coef[1,"Pr(>|t|)"], digits = 3))
  } else {
    SR.st[i,1] <- NA
    SR.st[i,2] <- NA
    SR.st[i,3] <- NA
  }
}

colnames(SR.st) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.st.2)[0,])
SR.st <- as.data.frame(cbind(congener, SR.st))

# Convert R2 and p-value to numeric
SR.st$`Sampling Rate (m3/d)` <- as.numeric(SR.st$`Sampling Rate (m3/d)`)
SR.st$R2 <- as.numeric(SR.st$R2)
SR.st$`p-value` <- as.numeric(SR.st$`p-value`)

# Update R2 and p-value to NA based on conditions
mask <- SR.st$R2 < 0.9 | SR.st$`p-value` > 0.05
SR.st$`Sampling Rate (m3/d)`[mask] <- NA
SR.st$R2[mask] <- NA
SR.st$`p-value`[mask] <- NA

# Values
SR.st.n <- length(na.omit(SR.st$`Sampling Rate (m3/d)`))
SR.st.ave <- mean(SR.st$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.st.sd <- sd(SR.st$`Sampling Rate (m3/d)`, na.rm = TRUE)

#export
write.csv(SR.st, file = "Output/Data/csv/WDSamplingRateStV2.csv")

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
  if (PUF.mean.2[i] > 0) {
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

# Values
SR.rot.n <- length(na.omit(SR.rot$`Sampling Rate (m3/d)`))
SR.rot.ave <- mean(SR.rot$`Sampling Rate (m3/d)`, na.rm = TRUE)
SR.rot.sd <- sd(SR.rot$`Sampling Rate (m3/d)`, na.rm = TRUE)

#export
write.csv(SR.rot, file = "Output/Data/csv/WDSamplingRateRotV2.csv")

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
  annotate("text", x = 1, y = 5.5, label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0, label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
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
ggsave("Output/Plots/SamplingRates/PCB18.30SamplingRates.png",
       plot = plot.18.30, width = 8, height = 10, dpi = 1300)

# PCB 52
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

# Plotting
plot.52 <- ggplot(combined_data, aes(x = time, y = PCB52_comb,
                                     color = group, fill = group)) +
  geom_point(shape = 21, size = 6, color = "black") +
  stat_smooth(method = "lm", se = FALSE, aes(group = group),
              formula = y ~ 0 + x, fullrange = TRUE) +
  annotate("text", x = 1, y = 5.5, label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0, label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  theme_bw() +
  xlim(0, 125) +
  ylim(0, 6) +
  xlab(expression(bold("Deployment time (h)"))) +
  ylab(expression(bold("Effective Volume PCB 52 (m"^"3"*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 24),
        axis.text.x = element_text(face = "bold", size = 22),
        axis.title.x = element_text(face = "bold", size = 24),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        aspect.ratio = 1.5)

# See plot
plot.52

# Save plot
ggsave("Output/Plots/SamplingRates/PCB52SamplingRates.png",
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
  annotate("text", x = 1, y = 5.5, label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0, label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
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
ggsave("Output/Plots/SamplingRates/PCB118SamplingRates.png",
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
  annotate("text", x = 1, y = 5.5, label = bquote("Dynamic" ~ "=" ~ .(round(slopes$slope[slopes$group == "Dynamic"], 2)) ~ "(m"^3*"/d)"),
           hjust = 0, size = 10, color = "black") +
  annotate("text", x = 1, y = 5.0, label = bquote("Static" ~ "=" ~ .(round(slopes$slope[slopes$group == "Static"], 2)) ~ "(m"^3*"/d)"),
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
ggsave("Output/Plots/SamplingRates/PCB187SamplingRates.png",
       plot = plot.187, width = 8, height = 10, dpi = 1300)
