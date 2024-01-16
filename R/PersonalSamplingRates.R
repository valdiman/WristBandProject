## Code to individual PCB calculate "personal" sampling rates
# for silicone wristbands. Sampling rates were calculated with

# Install packages
install.packages("readxl") #say no!
install.packages("gridExtra")
install.packages("ggplot2")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
}

# Read data from excel ----------------------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                             col_names = TRUE, col_types = NULL))
data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))

# Calculate personal sampling rate Amanda ---------------------------------
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
SR.amanda.r <- data.frame(SR.amanda.r, group = "amanda.r")

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
SR.amanda.l <- data.frame(SR.amanda.l, group = "amanda.l")

# Plot
# Convert numerical columns to numeric
SR.amanda.r[, 2:4] <- apply(SR.amanda.r[, 2:4], 2, as.numeric)
# Organize PCB names
SR.amanda.r$congener <- factor(SR.amanda.r$congener,
                            levels = unique(SR.amanda.r$congener))

# Define colors for each group
group_colors <- c("amanda.r" = "blue")

# Plot with legend
ggplot(SR.amanda.r[SR.amanda.r$Sampling_Rate > 0 & SR.amanda.r$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  scale_color_manual(values = group_colors, name = "Group") +
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

# Define colors for each group
group_colors <- c("amanda.l" = "red")

ggplot(SR.amanda.l[SR.amanda.l$Sampling_Rate > 0 & SR.amanda.l$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Calculate personal sampling rate Kay ------------------------------------
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
SR.kay.r <- data.frame(SR.kay.r, group = "kay.r")

# Plot
SR.kay.r[, 2:4] <- apply(SR.kay.r[, 2:4], 2, as.numeric)
# Organize PCB names
SR.kay.r$congener <- factor(SR.kay.r$congener,
                               levels = unique(SR.kay.r$congener))

# Define colors for each group
group_colors <- c("kay.r" = "green")

ggplot(SR.kay.r[SR.kay.r$Sampling_Rate > 0 & SR.kay.r$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates Right (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Combine sampling rates --------------------------------------------------
# Combine the three data frames
combined_SR <- rbind(SR.amanda.r, SR.amanda.l, SR.kay.r)

# Plot the combined data with different colors for each group
ggplot(combined_SR[combined_SR$Sampling_Rate > 0 & combined_SR$p_value < 0.05, ],
       aes(x = congener, y = Sampling_Rate, color = group)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Sampling Rates (m"^3*"/d)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 5, angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Calculate personal sampling rate Ya'u -----------------------------------
{
  # Select WBs to calculate air concentration
  data.yau.1 <- data.yau[1:6, 4:174]
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  time <- data.yau[1:6, 1]
  conc <- t(sweep(data.yau.1, 1, 0.5 * time, "/"))
  
  
  
  
  
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
SR.amanda.r <- data.frame(SR.amanda.r, group = "amanda.r")


