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
  bl <- read.csv("Data/WBTeBlanks.csv")
  wt <- read.csv("Data/WBteSamples.csv")
}

# Distribution analysis ---------------------------------------------------
# Remove metadata from blank data
bl.1 <- subset(bl, select = -c(sample.code))
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
  qqnorm(bl.1$PCB9, main = "Concentration (ng/g)")
  qqline(bl.1$PCB9)
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
# Create LOQ, i.e., upper 95 CI% (mean + 1.96*sd/(n)^0.5)
loq <- colMeans(bl.1) + 1.96*sapply(bl.1, sd)/sqrt(21)
loq <- data.frame(t(loq))

# Sample loq comparison ---------------------------------------------------
# If s.1 > loq, then wt[, 2:174], if not 0
# Remove sample names from wt
wt.1 <- wt[, 2:174]
# Create matrix to storage s.1 or loq values in s.2
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
# Add sample
wt.2 <- cbind(wt$code.teacher, wt.2)
# Change column name
names(wt.2)[names(wt.2) == 'wt$code.teacher'] <- 'code.teacher'
# Add column names
colnames(wt.2)[2:174] <- colnames(wt)[2:174]

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
  colnames(prof) <- wt$code.teacher
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
ggplot(prof, aes(x = congener, y = wt.01.r)) + # change y
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
costheta <- matrix(nrow = length(prof.cos[1,]),
                      ncol = length(prof.cos[1,]))

# Perform Cosine Theta
for (i in 1:length(prof.cos[1,])) {
  for (j in 1:length(prof.cos[1,])) {
    m1 <- prof.cos[,i]
    m2 <- prof.cos[,j]
    costheta[i,j] <- sum(m1*m2)/(sum(m1^2)*sum(m2^2))^0.5
  }
}

# Just 3 significant figures
costheta <- formatC(signif(costheta, digits = 3))
# Remove upper diagonal values
costheta[upper.tri(costheta)] <- NA
# Add name to columns
colnames(costheta) <- colnames(prof.cos)
# Add names to rows
rownames(costheta) <- colnames(prof.cos)
# Export data
write.csv(costheta, file = "Output/Data/csv/Teachers/costheta.csv")

# Predict Concentrations  ------------------------------------------------
# Read kos
ko.p <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSRV02.csv")
pcb_list <- ko.p$congener
# Read logKoa
logKoa <- read.csv("Data/logKoa.csv")
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
wt.3 <- cbind(wt$school.year, wt.3)
# Add time back
wt.3 <- cbind(wt$time.day, wt.3)
# Add sample
wt.3 <- cbind(wt$code.teacher, wt.3)
# Add WB volume
wt.3 <- cbind(wt$vol.WB, wt.3)
# Add WB area
wt.3 <- cbind(wt$area.WB, wt.3)
# Change column names
names(wt.3)[names(wt.3) == 'wt$school.year'] <- 'school.year'
names(wt.3)[names(wt.3) == 'wt$code.teacher'] <- 'code.teacher'
names(wt.3)[names(wt.3) == 'wt$time.day'] <- 'time.day'
names(wt.3)[names(wt.3) == 'wt$vol.WB'] <- 'vol.WB'
names(wt.3)[names(wt.3) == 'wt$area.WB'] <- 'area.WB'

# Calculate Veff
vol_matrix <- matrix(rep(wt.3$vol.WB, each = 171), nrow = 36, byrow = TRUE)
area_matrix <- matrix(rep(wt.3$area.WB, each = 171), nrow = 36, byrow = TRUE)
time_matrix <- matrix(rep(wt.3$time.day, each = 171), nrow = 36, byrow = TRUE)
logK_matrix <- matrix(rep(logKwb$logKwb, times = 36), nrow = 36, byrow = FALSE)
ko2 <- ko.p$Average_ko2

# Calculate veff.teacher as a 36 x 171 matrix
veff.teacher <- 10^logK_matrix * vol_matrix * 
  (1 - exp(-ko2 * area_matrix / vol_matrix / 10^logK_matrix * time_matrix))

# Estimate concentration from worn WBs
wt.mass <- wt.3[, 6:176]
conc.WB <- wt.mass / veff.teacher
conc.WB <- as.data.frame(conc.WB)
conc.WB$code.teacher <- wt$code.teacher
conc.WB$school.year <- wt$school.year
conc.WB$sample <- wt$Congener.Sample

# Export data
write.csv(conc.WB, file = "Output/Data/csv/Teachers/ConcentrationTeachers.csv")

# Predicted Total PCB Concentration ---------------------------------------
tPCB.conc.WB <- as.data.frame(rowSums(conc.WB[, 1:171], na.rm = TRUE))
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
  mutate(prefix = str_sub(code.teacher, 1, -3)) %>%
  group_by(prefix) %>%
  mutate(
    mean_tPCB = ifelse(n() > 1, mean(tPCB, na.rm = TRUE), tPCB),
    sd_tPCB = ifelse(n() > 1, sd(tPCB, na.rm = TRUE), NA),
    rsd_tPCB = ifelse(!is.na(sd_tPCB) & mean_tPCB != 0, (sd_tPCB / mean_tPCB) * 100, NA)
  ) %>%
  ungroup() %>%
  select(-prefix)  # Remove the temporary prefix column

# Remove duplicates and keep only one row for each prefix
plot_data <- updated_tPCB %>%
  group_by(prefix = str_sub(code.teacher, 1, -3)) %>%  # Group by prefix
  slice(1) %>%  # Keep only the first entry for each prefix group
  ungroup() %>%
  mutate(label = str_remove(prefix, "\\."))  # Remove dot and last character from label

# Plot the data
tPCB.wt <- ggplot(plot_data, aes(x = factor(label), y = mean_tPCB)) +
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  geom_errorbar(aes(ymin = mean_tPCB - sd_tPCB, ymax = mean_tPCB + sd_tPCB), width = 0.2) +
  theme_bw() +
  theme(aspect.ratio = 0.5) +
  ylab(expression(bold("Estimated Air Concentration " *Sigma*"PCB (ng/m3)"))) +
  xlab("") +
  theme(axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(face = "bold"))

# Display the plot
print(tPCB.wt)

# Save plot in folder
ggsave("Output/Plots/Teachers/WTtPCBVeff2.png", plot = tPCB.wt,
       width = 10, height = 5, dpi = 1200)

# Create a cleaned version of 'code.teacher' without the .l or .r
updated_tPCB <- updated_tPCB %>%
  mutate(code.teacher.clean = sub("\\.[lr]$", "", code.teacher))  # Remove the last character and dot

# Remove duplicate rows where the mean was calculated
unique_tPCB <- updated_tPCB %>%
  distinct(code.teacher.clean, .keep_all = TRUE)

# Plot the mean PCB concentrations against 'school.year'
ggplot(unique_tPCB, aes(x = school.year, y = mean_tPCB, label = code.teacher.clean)) +
  geom_point(size = 4, color = "black") +  # Plot the mean as points
  geom_text(vjust = -1, hjust = 0.5, fontface = "bold", size = 3.5) +  # Add cleaned 'code.teacher' labels
  theme_bw() +
  ylim(0, 10) +
  ylab(expression(bold("Estimated Air Concentration " *Sigma*"PCB (ng/m3)"))) +
  xlab("School Year") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 12)) +
  ggtitle("Mean PCB Concentration by School Year")

# Improve version
tPCB.wt.yr <- ggplot(unique_tPCB, aes(x = school.year, y = mean_tPCB)) +
  geom_point(size = 4, shape =21, color = "black", fill = "white") +  # Plot the mean as points
  theme_bw() +
  ylim(0, 10) +
  ylab(expression(bold("Estimated Air Concentration " * Sigma * "PCB (ng/m³)"))) +
  xlab("School Built") +
  theme(
    axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14)
  )

# Display the plot
print(tPCB.wt.yr)

# Save plot in folder
ggsave("Output/Plots/Teachers/WTtPCByrVeff2.png", plot = tPCB.wt.yr,
       width = 10, height = 5, dpi = 1200)

# Add year renovations to schools. This is too much speculation
# Change Lincoln Elementary from 1918 to 1974
unique_tPCB[2, 3] <- 1974
# Change Helen Lemme Elementary from 1970 to 2011
unique_tPCB[15, 3] <- 2011

# Plot the mean PCB concentrations against 'school.year'
ggplot(unique_tPCB, aes(x = school.year, y = mean_tPCB, label = code.teacher.clean)) +
  geom_point(size = 4, color = "black") +  # Plot the mean as points
  geom_text(vjust = -1, hjust = 0.5, fontface = "bold", size = 3.5) +  # Add cleaned 'code.teacher' labels
  theme_bw() +
  ylab(expression(Sigma*"PCB (ng/m3)")) +
  xlab("School Year") +
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 12))

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
plot.pcb.data <- conc.WB[, 1:171]

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
ggsave("Output/Plots/Teachers/WTPCBiVeff2.png", plot = plot.pcbi,
       width = 10, height = 5, dpi = 1200)

# Concentration Profile Analysis ------------------------------------------
# Profiles
tmp <- rowSums(conc.WB[, 1:171], na.rm = TRUE)
prof.WB.conc <- sweep(conc.WB[, 1:171], 1, tmp, FUN = "/")
prof.WB.conc <- cbind(conc.WB$code.teacher, prof.WB.conc)
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.WB.conc[, 2:172], na.rm = TRUE)

# See top PCBi for each sample
# Extract numeric columns and their names
numeric_data <- prof.WB.conc[, -1]
column_names <- colnames(numeric_data)

# Function to get the names of the top 10 columns
top_10_PCBi <- function(row) {
  top_indices <- order(row, decreasing = TRUE)[1:10]
  column_names[top_indices]
}

# Apply the function to each row
top_10_PCBi_matrix <- t(apply(numeric_data, 1, top_10_PCBi))

# Convert to dataframe and add the non-numeric column
top_10_PCBi_df <- as.data.frame(top_10_PCBi_matrix)
top_10_PCBi_df <- cbind(prof.WB.conc[, 1], top_10_PCBi_df)

# Rename columns
colnames(top_10_PCBi_df) <- c("code.teacher", paste0("Top10_PCBi_", 1:10))

# Add school built year
top_10_PCBi_df <- cbind(conc.WB$school.year, top_10_PCBi_df)

# Create the plot
ggplot(top_10_PCBi_df, aes(x = Top10_PCBi_1, y = `conc.WB$school.year`)) +
  geom_point(color = "blue") +  # Scatter plot with blue points
  labs(x = "School Year", y = "Top 10 PCB Concentration (Top10_PCBi_1)",
       title = "School Year vs. Top 10 PCB Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# PCA
t.prof.WB.conc <- data.frame(t(prof.WB.conc[, 2:172]))
colnames(t.prof.WB.conc) <- conc.WB$code.teacher
PCA <- prcomp(t.prof.WB.conc, scale. = TRUE)
summary(PCA)

# Plot PCA results
ggplot2::autoplot(PCA, data = t.prof.WB.conc, size = 1.5, shape = 21) +
  theme_bw() +
  xlim(-0.4, 0.6) +
  ylim(-0.4, 0.6) +
  theme(aspect.ratio = 10/10)

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
write.csv(costheta_df, file = "Output/Data/csv/Teachers/CosineThetaTeachersVeff2.csv")

# Cosine theta visualization ----------------------------------------------
# Calculate tPCB values to be added to the plot
tPCB.conc.WB <- as.data.frame(rowSums(conc.WB[, 1:171], na.rm = TRUE))
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

# Plot the data
plot.cos.theta.low <- ggplot(data = costheta_correlations[1:13, ], # low values =< 0.14
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
ggsave("Output/Plots/Profiles/Teachers/CosThetaLowVeff2.png", plot = plot.cos.theta.low,
       width = 10, height = 10, dpi = 1200)

plot.cos.theta.high <- ggplot(data = costheta_correlations[586:595, ], # high values >= 0.86
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
ggsave("Output/Plots/Profiles/Teachers/CosThetaHighVeff2.png", plot = plot.cos.theta.high,
       width = 10, height = 10, dpi = 1200)

# Plot Individual PCB Profiles --------------------------------------------
# Select rows for the sample
selected_rows <- prof.WB.conc %>%
  #wt.19.l = Teacher 1
  #wt.25.r = Teacher 2
  filter(conc.WB$code.teacher %in% c("wt.25.r")) # need to change the sample!

# Create Source vector with correct length and values
source_vector <- rep(selected_rows[[1]], each = ncol(selected_rows) - 1)

# Reshape data from wide to long format, excluding non-PCB columns
prof_combined <- selected_rows %>%
  pivot_longer(cols = starts_with("PCB"), names_to = "congener", values_to = "Conc") %>%
  mutate(congener = sub("PCB", "", congener))

# Add Source to the reshaped data
prof_combined <- prof_combined %>%
  mutate(Source = source_vector)

prof_combined$congener <- gsub("\\.", "+", prof_combined$congener)

# Convert congener to factor and set levels
prof_combined$congener <- factor(prof_combined$congener, 
                                 levels = unique(prof_combined$congener))

# Create the bar plot
plot.25.r <- ggplot(prof_combined, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 0.9, 
           color = "black", linewidth = 0.2) +
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(
    values = c("wt.25.r" = "blue"), 
    labels = c("wt.25.r" = "Teacher 2"),
    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(e)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

plot.16.r
plot.19.l
plot.25.r

# Save plot
ggsave("Output/Plots/Profiles/Teachers/Teacher2.png", plot = plot.25.r,
       width = 10, height = 3, dpi = 1200)

