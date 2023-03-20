## Code to individual PCB calculate sampling rates
# for silicone wristbands. Sampling rates were calculated with
# (i) static (2 times) and (ii) rotating set up.

# Install packages
install.packages("readxl") #say no!
install.packages("gridExtra")
install.packages("ggplot2")

# Load libraries
library(readxl)
library(ggplot2)
library(gridExtra)

# Sampling rates calculations under static conditions ---------------------

# Read data.xlsx
PUF <- data.frame(read_excel("DataStatic.xlsx", sheet = "PUFv2",
                               col_names = TRUE, col_types = NULL))

WB <- data.frame(read_excel("DataStatic.xlsx", sheet = "WDv2",
                                col_names = TRUE, col_types = NULL))

# Remove metadata
WB.1 <- subset(WB, select = -c(sample:time))

# Perform linear regression of WD accumulate mass vs time to check linearity
# create matrix
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
write.csv(WBMatrix, file = "WDLinearity.csv")

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
    SR[i,1] <- 0
    SR[i,2] <- 0
    SR[i,3] <- 0
  }
}

colnames(SR) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.1)[0,])
SR <- cbind(congener, SR)

#export
write.csv(SR, file = "WDSamplingRate.csv")

# Plots
# Change number congener in [] 

fit1 <- lm(WB.1[, 152]/PUF.mean[152] ~ 0 + time)
ggplot(WB, aes(x = time, y = WB.1[, 152]/PUF.mean[152])) +
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
  theme(plot.title = element_text(size = 8, face = "bold")) +
  annotate("text", x = 25, y = 3,
           label = paste(" Slope (m3/d)=",
                         signif(summary(fit1)$coef[1,"Estimate"], 2),
                         ", R2 = ", signif(summary(fit1)$adj.r.squared, 3)),
                         size = 2.5, fontface = 2)

# Sampling rates calculations under static and rotating conditions -------------------

# Read data.xlsx
PUF.v2 <- data.frame(read_excel("DataV04.xlsx", sheet = "PUFv3",
                             col_names = TRUE, col_types = NULL))

WB.st <- data.frame(read_excel("DataV03.xlsx", sheet = "WDST",
                            col_names = TRUE, col_types = NULL))

WB.rot <- data.frame(read_excel("DataV03.xlsx", sheet = "WDROT",
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
    SR.st[i,1] <- 0
    SR.st[i,2] <- 0
    SR.st[i,3] <- 0
  }
}

colnames(SR.st) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.st.2)[0,])
SR.st <- cbind(congener, SR.st)

#export
write.csv(SR.st, file = "WDSamplingRateStV3.csv")

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
  annotate("text", x = 25, y = 3,
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
    SR.rot[i,1] <- 0
    SR.rot[i,2] <- 0
    SR.rot[i,3] <- 0
  }
}

colnames(SR.rot) <-c("Sampling Rate (m3/d)", "R2", "p-value")
congener <- names(head(WB.rot.2)[0,])
SR.rot <- cbind(congener, SR.rot)

#export
write.csv(SR.rot, file = "WDSamplingRateRotV2.csv")

# Plots
# Change number congener in [] 

fit1 <- lm(WB.rot.2[, 153]/PUF.mean.2[153] ~ 0 + time.2)
ggplot(WB.rot.2, aes(x = time.2*24, y = WB.rot.2[, 153]/PUF.mean.2[153])) +
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
  ylim(0, 6) +
  theme(plot.title = element_text(size = 8, face = "bold")) +
  annotate("text", x = 25, y = 6,
           label = paste(" Slope (m3/d)=",
                         signif(summary(fit1)$coef[1,"Estimate"], 2),
                         ", R2 = ", signif(summary(fit1)$adj.r.squared, 3)),
           size = 2.5, fontface = 2)

