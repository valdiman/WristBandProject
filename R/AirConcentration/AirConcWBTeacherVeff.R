## Script to analyze school teacher wristband data
# estimate concentrations and PCB profiles

# Install packages
install.packages("reshape2")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("tibble")
install.packages("ggfortify")
install.packages("stringr")

# Library
{
  library(ggplot2)
  library(reshape2) # For melt function
  library(tidyr) # Data manipulation
  library(dplyr) # performs %>%
  library(tibble) # add column
  library(ggfortify) # PCA analysis
  library(stringr)
}

# Read data ---------------------------------------------------------------
{
  bl.0 <- read.csv("Data/IRO/BlankMassStudy3_4_5.csv")
  bl <- bl.0[8:28, c(1, 8:180)] # select rows and columns
  wt.0 <- read.csv("Data/IRO/SampleMassStudy3_4_5.csv")
  wt <- wt.0[39:74, c(1, 7, 9:10, 12:184)]
}

# Distribution analysis ---------------------------------------------------
# Remove metadata from blank data
bl.1 <- subset(bl, select = -c(sid))
# Look at the distribution of the blank data
# Create matrix to storage data
normality <- matrix(nrow = length(bl.1[1,]), ncol = 2)

for (i in 1:length(bl.1[1, ])) {
  normality[i, 1] <- shapiro.test(bl.1[,i])$p.value
  normality[i, 2] <- shapiro.test(log10(bl.1[,i]))$p.value
}

# Just 3 significant figures
normality <- formatC(signif(normality, digits = 3))
# Add congener names
congeners <- colnames(bl.1)
# Include congener names
normality <- cbind(congeners, normality)
# Change column names
colnames(normality) <- c("Congener", "shapiro.normal", "shapiro.log10")
# Create Q-Q plot for individual PCB congeners
{
  qqnorm(bl.1$PCB52, main = "Concentration (ng/g)")
  qqline(bl.1$PCB52)
}

# Bar plot to visualize Shapiro test (normality)
# Organize data to be plotted
norma <- data.frame(normality)
norma$Congener <- as.character(norma$Congener)
norma$shapiro.normal<- as.numeric(as.character(norma$shapiro.normal))
norma$shapiro.log10 <- as.numeric(as.character(norma$shapiro.log10))

# Plot
# Values less than -log10(0.05) are not significant
# Normal data
ggplot(norma, aes(x = Congener, y = -log10(shapiro.normal))) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "black") +
  theme_test() +
  theme(aspect.ratio = 2/12) +
  ylab("-log10 (p-value)") +
  theme(axis.title.y = element_text(face = "bold", size = 8)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "red",
             linewidth = 0.8)

# log10 data
ggplot(norma, aes(x = Congener, y = -log10(shapiro.log10))) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "black") +
  theme_test() +
  theme(aspect.ratio = 2/12) +
  ylab("-log10 (p-value)") +
  theme(axis.title.y = element_text(face = "bold", size = 8)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "red",
             linewidth = 0.8)
  
# Export data
write.csv(normality, file = "Output/Data/csv/Teachers/BlankNormalityTeachers.csv")

# Calculate LOQ -----------------------------------------------------------
# From above analysis, normal scale is used to calculate LOQ
# Create LOQ, i.e., upper 95 CI% (mean + 1.96*sd/(n)^0.5)
loq <- colMeans(bl.1) + 1.96*sapply(bl.1, sd)/sqrt(21)
loq <- data.frame(t(loq))

# Sample loq comparison ---------------------------------------------------
# If s.1 > loq, then wt[, 5:177], if not 0
# Remove metadata from wt
wt.1 <- wt[, 5:177]
# Create matrix to storage
wt.2 <- matrix(NA, nrow = dim(wt.1)[1], ncol = dim(wt.1)[2])
# Do comparison
for(i in 1:dim(wt.1)[1]) {
  for(j in 1:dim(wt.1)[2]) {
    if (wt.1[i, j] > loq[j]) {
      wt.2[i, j] <- wt.1[i, j]
    } else {
      wt.2[i, j] <- 0
    }
  }
}

# Transform to data.frame
wt.2 <- data.frame(wt.2)
# Add sample code
wt.2 <- cbind(wt$sid, wt.2)
# Change column name
names(wt.2)[names(wt.2) == 'wt$sid'] <- 'code.teacher'
# Add column names (congener)
colnames(wt.2)[2:174] <- colnames(wt)[5:177]

# Individual PCB detection frequency --------------------------------------
wt.fr <- as.data.frame(colSums(wt.2[2:174] > 0)/length(wt.2[, 1])*100)
colnames(wt.fr) <- c("Detection.Frequency")

wt.fr.stats <- wt.fr %>%
  summarize(
    mean_detection = mean(Detection.Frequency),
    sd_detection = sd(Detection.Frequency),
    min_detection = min(Detection.Frequency),
    max_detection = max(Detection.Frequency)
  )

# Add PCB names as a new column in wt.fr
wt.fr$congener <- colnames(wt.2)[2:174]
wt.fr$congener <- gsub('\\.', '+', wt.fr$congener) # replace dot for +
wt.fr$congener <- factor(wt.fr$congener,
                                 levels = rev(wt.fr$congener)) # change the order to be plotted.

# Detection frequency plot
plot.cong.freq <- ggplot(wt.fr, aes(x = Detection.Frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "white", color = "black") +
  ylab("") +
  theme_bw() +
  xlim(c(0, 100)) +  # Set x-axis limits to 0-100%
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(aspect.ratio = 20/5) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.text.y = element_text(face = "bold", size = 6))

# Display the plot
print(plot.cong.freq)

# Save plot in folder
ggsave("Output/Plots/Teachers/FreqPCBWT.png", plot = plot.cong.freq,
       width = 5, height = 10, dpi = 1200)

# Plot tPCB ---------------------------------------------------------------
ggplot(wt.2, aes(y = rowSums(wt.2[2:174]), x = factor(code.teacher))) +
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression(Sigma*"PCB (ng/g)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1, vjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Mass Profile Analysis --------------------------------------------------------
# Generate profile
{
  tmp <- rowSums(wt.1, na.rm = TRUE)
  prof <- sweep(wt.1, 1, tmp, FUN = "/")
  # Transpose
  prof <- t(prof)
  # Include sample names
  colnames(prof) <- wt$sid
  # Get congeners
  congener <- row.names(prof)
  # Add congener names to first column
  prof <- cbind(congener, prof)
  # Delete row names
  rownames(prof) <- NULL
  # Transform prof column to data.frame
  prof <- data.frame(prof)
  # Change "." to "+" in the PCB congeners
  prof$congener <- lapply(prof$congener, function(x) {gsub("\\.", "+", x)})
  # Convert all character columns to numeric
  prof[, 2:37] <- as.data.frame(apply(prof[, 2:37], 2, as.numeric))
  # Organize PCB congeners
  prof$congener <- factor(prof$congener, levels = unique(prof$congener)) 
}

# Plots
ggplot(prof, aes(x = congener, y = S02_1)) + # change y (sid)
  geom_bar(position = position_dodge(), stat = "identity",
           fill = "black") +
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 4/16) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Remove congener column
prof.cos <- prof[, 2:37]
# Create matrix to storage results
costheta.mass <- matrix(nrow = length(prof.cos[1,]),
                      ncol = length(prof.cos[1,]))

# Perform Cosine Theta
for (i in 1:length(prof.cos[1,])) {
  for (j in 1:length(prof.cos[1,])) {
    m1 <- prof.cos[,i]
    m2 <- prof.cos[,j]
    costheta.mass[i,j] <- sum(m1*m2)/(sum(m1^2)*sum(m2^2))^0.5
  }
}

# Just 3 significant figures
costheta.mass <- formatC(signif(costheta.mass, digits = 3))
# Remove upper diagonal values
costheta.mass[upper.tri(costheta.mass)] <- NA
# Add name to columns
colnames(costheta.mass) <- colnames(prof.cos)
# Add names to rows
rownames(costheta.mass) <- colnames(prof.cos)
# Export data
write.csv(costheta.mass, file = "Output/Data/csv/Teachers/CosineThetaMass.csv")

# Predict Concentrations  ------------------------------------------------
# Read kos
ko.p <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSR.csv")
pcb_list <- ko.p$congener
# Read logKoa
logKoa <- read.csv("Data/IRO/logKoa.csv")
# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# Subset logKwb to include same congeners as ko
logKwb <- logKwb[logKwb$congener %in% pcb_list, ]
# Get ko
ko.p <- ko.p[7]
# Subset wt.1
wt.3 <- wt.1[, intersect(pcb_list, colnames(wt.1))]
# Add school year
wt.3$school.year <- c(1968, 1968, 1926, 1926, NA, 2017, 1962, 1939,
  1939, 1918, 1918, 1939, 1939, 1939, 1939, 2017,
  2006, 1972, 1972, 1970, 1970, 1968, 1968, 1968,
  1968, 1968, 1968, 1968, 2017, 2017, 2017, 2017,
  2017, 2017, 2017, 2017)

# Add time back
wt.3 <- cbind(wt$time / 24, wt.3) # transform hours to day
# Add sample
wt.3 <- cbind(wt$sid, wt.3)
# Add WB volume
wt.3 <- cbind(wt$vol.WB, wt.3)
# Add WB area
wt.3 <- cbind(wt$area.WB, wt.3)
# Change column names
names(wt.3)[names(wt.3) == 'wt$sid'] <- 'code.teacher'
names(wt.3)[names(wt.3) == 'wt$time/24'] <- 'time.day'
names(wt.3)[names(wt.3) == 'wt$vol.WB'] <- 'vol.WB'
names(wt.3)[names(wt.3) == 'wt$area.WB'] <- 'area.WB'

n_teachers <- nrow(wt.3)
n_pcbs <- nrow(logKwb)

# Calculate Veff
vol_matrix  <- matrix(rep(wt.3$vol.WB, times = n_pcbs),
                      nrow = n_teachers, ncol = n_pcbs)
area_matrix <- matrix(rep(wt.3$area.WB, times = n_pcbs),
                      nrow = n_teachers, ncol = n_pcbs)
time_matrix <- matrix(rep(wt.3$time.day, times = n_pcbs),
                      nrow = n_teachers, ncol = n_pcbs)
logK_matrix <- matrix(rep(logKwb$logKwb, each = n_teachers),
                      nrow = n_teachers, ncol = n_pcbs)
ko2_matrix  <- matrix(rep(ko2, each = n_teachers),
                      nrow = n_teachers, ncol = n_pcbs)

# Calculate veff.teacher as a 36 x 173 matrix
veff.teacher <- 10^logK_matrix * vol_matrix *
  (1 - exp(-ko2_matrix * area_matrix / vol_matrix / 10^logK_matrix * time_matrix))

# Estimate concentration from worn WBs
wt.mass <- wt.3[, 5:177]
conc.WB <- wt.mass / veff.teacher
conc.WB <- as.data.frame(conc.WB)
conc.WB$code.teacher <- wt.3$code.teacher
conc.WB$school.year <- wt.3$school.year

# Export data
write.csv(conc.WB,
          file = "Output/Data/csv/Teachers/ConcentrationTeachers.csv",
          row.names = FALSE)

# Predicted Total PCB Concentration ---------------------------------------
tPCB.conc.WB <- as.data.frame(rowSums(conc.WB[, 1:173], na.rm = TRUE))
tPCB.conc.WB$code.teacher <- conc.WB$code.teacher
tPCB.conc.WB$school.year <- conc.WB$school.year
# Change columns names
colnames(tPCB.conc.WB) <- c("tPCB", "code.teacher", "school.year")

# See
tPCB.conc.WB
min(tPCB.conc.WB$tPCB)
max(tPCB.conc.WB$tPCB)
mean(tPCB.conc.WB$tPCB)
sd(tPCB.conc.WB$tPCB)

# Need to change format to include duplicate
tPCB.conc.WB <- tPCB.conc.WB %>%
  arrange(code.teacher) %>%
  mutate(ID = row_number())

ggplot(tPCB.conc.WB, aes(x = factor(ID), y = tPCB, fill = code.teacher)) +
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_x_discrete(labels = tPCB.conc.WB$code.teacher) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression(Sigma*"PCB (ng/m3)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7))

# Averaging duplicates
updated_tPCB <- tPCB.conc.WB %>%
  mutate(prefix = str_remove(code.teacher, "_[^_]+$")) %>%
  group_by(prefix) %>%
  mutate(
    mean_tPCB = ifelse(n() > 1, mean(tPCB, na.rm = TRUE), tPCB),
    sd_tPCB   = ifelse(n() > 1, sd(tPCB, na.rm = TRUE), NA),
    rsd_tPCB  = ifelse(!is.na(sd_tPCB) & mean_tPCB != 0,
                       (sd_tPCB / mean_tPCB) * 100, NA)
  ) %>%
  ungroup() %>%
  select(-prefix)

# Remove duplicates and keep only one row for each prefix
plot_data <- updated_tPCB %>%
  mutate(prefix = str_remove(code.teacher, "_[^_]+$")) %>%
  group_by(prefix) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    label = prefix)

# Plot the data
tPCB.wt <- ggplot(plot_data, aes(x = factor(label), y = mean_tPCB)) +
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  geom_errorbar(aes(ymin = mean_tPCB - sd_tPCB, ymax = mean_tPCB + sd_tPCB),
                width = 0.2) +
  theme_bw() +
  theme(aspect.ratio = 0.5) +
  ylab(expression(bold("Estimated Air Concentration " *Sigma*"PCB (ng/m"^3 *")"))) +
  xlab("") +
  theme(axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold"))

# Display the plot
print(tPCB.wt)

# Save plot in folder
ggsave("Output/Plots/Teachers/WTtPCBVeff.png", plot = tPCB.wt,
       width = 10, height = 5, dpi = 1200)

# Plot the mean PCB concentrations against 'school.year'
ggplot(updated_tPCB, aes(x = school.year, y = mean_tPCB, label = code.teacher)) +
  geom_point(size = 4, color = "black") +  # Plot the mean as points
  geom_text(vjust = -1, hjust = 0.5, fontface = "bold", size = 3.5) +
  theme_bw() +
  ylim(0, 10) +
  ylab(expression(bold("Estimated Air Concentration " *Sigma*"PCB (ng/m3)"))) +
  xlab("School Year") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 12)) +
  ggtitle("Mean PCB Concentration by School Year")

# Improve version
tPCB.wt.yr <- ggplot(updated_tPCB, aes(x = school.year, y = tPCB)) +
  geom_point(size = 4, shape = 21, color = "black", fill = "white") +
  theme_bw() +
  ylab(expression(bold("Estimated Air Concentration " * Sigma * "PCB (ng/m³)"))) +
  xlab("School Built") +
  theme(
    axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14))

# Display the plot
print(tPCB.wt.yr)

# Save plot in folder
ggsave("Output/Plots/Teachers/WTtPCByrVeff.png", plot = tPCB.wt.yr,
       width = 10, height = 5, dpi = 1200)

# Predicted PCBi Concentrations -------------------------------------------
# Calculate means and standard deviations for columns 1 to 173
median_conc <- apply(conc.WB[, 1:171], 2, median, na.rm = TRUE)
mean_conc <- colMeans(conc.WB[, 1:171], na.rm = TRUE)
st_conc <- apply(conc.WB[, 1:171], 2, sd, na.rm = TRUE)
cov_conc <- st_conc / mean_conc * 100

PCB.conc.WB.ave <- data.frame(
  Column = names(mean_conc),
  Mean_Conc = mean_conc,
  St_Conc = st_conc,
  COV = cov_conc,
  Median_Conc = median_conc
)

min(PCB.conc.WB.ave$Mean_Conc)
max(PCB.conc.WB.ave$Mean_Conc)
mean(PCB.conc.WB.ave$COV)

# Top 10
top_10_PCB.means <- PCB.conc.WB.ave %>%
  arrange(desc(Mean_Conc)) %>%
  head(10)
print(top_10_PCB.means)

# Plot PCBi
plot.pcb.data <- conc.WB[, 1:173]

# Reshape the data from wide to long format
plot.pcb.long <- pivot_longer(
  plot.pcb.data, 
  cols = everything(),          # Convert all congener columns
  names_to = "Congener",         # Column name for congener names
  values_to = "Concentration"    # Column name for concentration values
)

plot.pcb.long <- plot.pcb.long %>%
  mutate(Congener = factor(Congener, levels = colnames(plot.pcb.data)))

plot.pcbi <- ggplot(plot.pcb.long, aes(x = Congener, y = Concentration)) +
  geom_boxplot(outlier.shape = 21, outlier.colour = "black", outlier.fill = NA) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 3/9) +
  ylab(expression(bold("Estimated Air PCBi Concentration (ng/m³)"))) +
  theme(
    axis.text.x = element_text(face = "bold", size = 5, angle = 90, hjust = 1, vjust = 0.6),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12)
  )

# Display the plot
plot.pcbi

# Save plot in folder
ggsave("Output/Plots/Teachers/WTPCBiVeff.png", plot = plot.pcbi,
       width = 10, height = 5, dpi = 1200)

# Concentration Profile Analysis ------------------------------------------
# Profiles
tmp <- rowSums(conc.WB[, 1:173], na.rm = TRUE)
prof.WB.conc <- sweep(conc.WB[, 1:173], 1, tmp, FUN = "/")
prof.WB.conc <- cbind(conc.WB$code.teacher, prof.WB.conc)
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.WB.conc[, 2:174], na.rm = TRUE)
t.prof.WB.conc <- data.frame(t(prof.WB.conc[, 2:174]))
colnames(t.prof.WB.conc) <- conc.WB$code.teacher

# Cosine theta analysis
# Create matrix to storage results
costheta.samples <- matrix(nrow = length(t.prof.WB.conc[1, ]),
                           ncol = length(t.prof.WB.conc[1, ]))

for (i in 1:length(t.prof.WB.conc[1,])) {
  for (j in 1:length(t.prof.WB.conc[1,])) {
    m1 <- t.prof.WB.conc[, i]
    m2 <- t.prof.WB.conc[, j]
    costheta.samples[i, j] <- sum(m1*m2)/(sum(m1^2)*sum(m2^2))^0.5
  }
}

# Remove upper diagonal values
costheta.samples[upper.tri(costheta.samples)] <- NA
# Add column names with samples name
colnames(costheta.samples) <- colnames(t.prof.WB.conc)
# Add row names with samples name
rownames(costheta.samples) <- colnames(t.prof.WB.conc)
costheta_df <- as.data.frame(costheta.samples)
# Add conc.wt$school.year as a new column to the data frame
costheta_df$school.year <- conc.WB$school.year

# Export
write.csv(costheta_df,
          file = "Output/Data/csv/Teachers/CosineThetaTeachersVeff.csv")

# Cosine theta visualization ----------------------------------------------
# Calculate tPCB values to be added to the plot
tPCB.conc.WB <- as.data.frame(rowSums(conc.WB[, 1:173], na.rm = TRUE))
colnames(tPCB.conc.WB) <- "tPCB"

# Add the tPCB values to costheta_df
costheta_df <- bind_cols(costheta_df, tPCB.conc.WB)

# Ensure column names are unique
colnames(costheta_df) <- make.unique(colnames(costheta_df))

# Extract the years and tPCB values
years <- costheta_df$school.year
tPCB <- costheta_df$tPCB

# Reshape the data from wide to long format, excluding 'school.year' and 'tPCB'
costheta_long <- costheta_df %>%
  select(-school.year, -tPCB) %>%
  pivot_longer(cols = everything(), names_to = "Var1", values_to = "Value") %>%
  mutate(Var2 = rep(colnames(costheta_df %>% select(-school.year, -tPCB)), each = nrow(costheta_df)))

# Define the reference column names (excluding year and tPCB)
ref_names <- colnames(costheta_df)[-c(37, 38)]

# Safely map year and tPCB based on Var1 and Var2 matching ref_names
costheta_long <- costheta_long %>%
  mutate(
    Var1_year = years[match(Var1, ref_names)],
    Var2_year = years[match(Var2, ref_names)],
    Var1_tPCB = round(tPCB[match(Var1, ref_names)], 2),
    Var2_tPCB = round(tPCB[match(Var2, ref_names)], 2),
    Var1 = paste0(Var1, " (", Var1_year, ", ", Var1_tPCB, ")"),
    Var2 = paste0(Var2, " (", Var2_year, ", ", Var2_tPCB, ")")
  ) %>%
  select(-Var1_year, -Var2_year, -Var1_tPCB, -Var2_tPCB)

# Remove NAs text from the long data
costheta_long <- costheta_long %>%
  filter(
    !str_detect(Var1, "NA"),
    !str_detect(Value, "NA"),
    !str_detect(Var2, "NA")
  )

# Sort the long data by correlation value in ascending order
costheta_correlations <- costheta_long %>%
  arrange(Value) %>%
  filter(Var1 != Var2)  # Exclude self-correlations (diagonal values)

# Change wt to tea
costheta_correlations <- costheta_correlations %>%
  mutate(
    Var1 = str_replace(Var1, "^wt", "tea"),
    Var2 = str_replace(Var2, "^wt", "tea")
  )

# Plot the data
plot.cos.theta.low <- ggplot(data = costheta_correlations[1:16, ], # low values =< 0.2
                             aes(x = Var1, y = Var2, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(face = "bold", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10)) +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 3) 

# Display the plot
plot.cos.theta.low

# Save plot
ggsave("Output/Plots/Profiles/Teachers/CosThetaLowVeff.png",
       plot = plot.cos.theta.low, width = 10, height = 10, dpi = 1200)

plot.cos.theta.high <- ggplot(data = costheta_correlations[583:595, ], # high values >= 0.85
                              aes(x = Var1, y = Var2, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(face = "bold", size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10)) +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 3)

# Display the plot
plot.cos.theta.high

# Save plot
ggsave("Output/Plots/Profiles/Teachers/CosThetaHighVeff.png",
       plot = plot.cos.theta.high, width = 10, height = 10, dpi = 1200)

# Plot Individual PCB Profiles --------------------------------------------
hand_vec <- wt.0$hand[39:74]

prof_long <- prof.WB.conc %>%
  # First, create a standard column "Source" with hand info
  mutate(Source = `conc.WB$code.teacher`,
         Source = sub("_\\d$", "", Source),
         Source = paste0(Source, ".", hand_vec),
         Source = paste0("Teacher ", Source)) %>%
  pivot_longer(
    cols = starts_with("PCB"),
    names_to = "congener",
    values_to = "Conc"
  ) %>%
  mutate(
    congener = gsub("\\.", "+", congener),
    congener = factor(congener, levels = unique(congener)))

# 2 PCB profiles for main text
prof_long_sel <- prof_long %>%
  filter(Source %in% c("Teacher S19.l", "Teacher S25.l")) %>%
  mutate(Source_label = Source)

plot.tea <- ggplot(prof_long_sel, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1,
           color = "black", linewidth = 0.2) +
  facet_wrap(~ Source_label, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.2), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " * Sigma * "PCB"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 9, face = "bold"
    ),
    axis.ticks.x.bottom = element_line())

plot.tea

# Save plot
ggsave("Output/Plots/Profiles/Teachers/Teacher2.png", plot = plot.tea,
       width = 22, height = 5, dpi = 500)

# All PCB profiles for SI
# "Teacher S01.l", "Teacher S01.r", "Teacher S02.l", "Teacher S02.r",
# "Teacher S03.r", "Teacher S04.r", "Teacher S05.r", "Teacher S06.l",
# "Teacher S06.r", "Teacher S07.l", "Teacher S07.r", "Teacher S08.r",
# "Teacher S09.r", "Teacher S10.r", "Teacher S11.r", "Teacher S12.r",
# "Teacher S13.r", "Teacher S14.r", "Teacher S15.r", "Teacher S15.l",
# "Teacher S16.r", "Teacher S17.r", "Teacher S18.l", "Teacher S18.r",
# "Teacher S20.r", "Teacher S21.r", "Teacher S22.r", "Teacher S23.l",
# "Teacher S23.r", "Teacher S24.l", "Teacher S24.r", "Teacher S25.r"

prof_long_sel <- prof_long %>%
  filter(Source %in% c("Teacher S23.r", "Teacher S24.l", "Teacher S24.r", "Teacher S25.r")) %>%
  mutate(Source_label = Source)  # optional, for faceting

plot.tea <- ggplot(prof_long_sel, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source_label, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.2), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " * Sigma * "PCB"))) +
  theme(
    axis.text.y = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 24, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 9, face = "bold"
    ),
    axis.ticks.x.bottom = element_line())

plot.tea

# Save plot
ggsave("Output/Plots/Profiles/Teachers/Teacher29-32.png", plot = plot.tea,
       width = 22, height = 15, dpi = 500)

