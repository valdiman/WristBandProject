## Concentration estimation for Study 3
# Office only

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("stringr")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(stringr)
}

# Read data ---------------------------------------------------------------
{
  data.s3 <- read.csv("Data/IRO/09_SampleWBMassStudy3.csv", check.names = FALSE)
  # Read individual PCB logKoa
  logKoa <- read.csv("Data/IRO/14_logKoa.csv")
  # ko from SamplingRates_ko.R file
  ko <- read.csv("Output/Data/Studies1_2/Study1/WDSamplingRateStatV1.csv")
  # Select only ko [m/d]
  ko <- ko[c(2,6)]
  # Modify "." to "+"
  ko$congener <- gsub("\\.", "+", ko$congener)
}

# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# Select only congeners and check that all have the same congeners
common_ids <- intersect(ko$congener, logKwb$congener)
ko <- ko[ko$congener %in% common_ids, ]
logKwb <- logKwb[logKwb$congener %in% common_ids, ]
data.s3.1 <- data.s3[, common_ids]

# Volunteers 4 to 14 (V)
# Calculate air PCB concentration from static WBs -------------------------
{
  # Calculate effective volume for static WBs
  # Use effective volume. Adult WBs
  Vwb.V4.V5 <- data.s3$vol.WB[1] # [m3]
  Awb.V4.V5 <- data.s3$area.WB[1] # [m2]
  veff_static.V4 <- 10^(logKwb$logKwb) * Vwb.V4.V5 * 
    (1 - exp(-ko$ko * Awb.V4.V5 / Vwb.V4.V5 / 10^(logKwb$logKwb) * data.s3$time[1] / 24))
  veff_static.V5 <- 10^(logKwb$logKwb) * Vwb.V4.V5 * 
    (1 - exp(-ko$ko * Awb.V4.V5 / Vwb.V4.V5 / 10^(logKwb$logKwb) * data.s3$time[1] / 24))
  Vwb.V6 <- data.s3$vol.WB[2] # [m3]
  Awb.V6 <- data.s3$area.WB[2] # [m2]
  veff_static.V6 <- 10^(logKwb$logKwb) * Vwb.V6 * 
    (1 - exp(-ko$ko * Awb.V6 / Vwb.V6 / 10^(logKwb$logKwb) * data.s3$time[2] / 24))
  # 3 replicated for volunteer 7
  Vwb.V7 <- data.s3$vol.WB[3] # [m3]
  Awb.V7 <- data.s3$area.WB[3] # [m2]
  veff_static.V7.1 <- 10^(logKwb$logKwb) * Vwb.V7 * 
    (1 - exp(-ko$ko * Awb.V7 / Vwb.V7 / 10^(logKwb$logKwb) * data.s3$time[3] / 24))
  veff_static.V7.2 <- 10^(logKwb$logKwb) * Vwb.V7 * 
    (1 - exp(-ko$ko * Awb.V7 / Vwb.V7 / 10^(logKwb$logKwb) * data.s3$time[4] / 24))
  veff_static.V7.3 <- 10^(logKwb$logKwb) * Vwb.V7 * 
    (1 - exp(-ko$ko * Awb.V7 / Vwb.V7 / 10^(logKwb$logKwb) * data.s3$time[5] / 24))
  Vwb.V8 <- data.s3$vol.WB[6] # [m3]
  Awb.V8 <- data.s3$area.WB[6] # [m2]
  veff_static.V8 <- 10^(logKwb$logKwb) * Vwb.V8 * 
    (1 - exp(-ko$ko * Awb.V8 / Vwb.V8 / 10^(logKwb$logKwb) * data.s3$time[6] /24))
  # 3 replicates for volunteer 9, but the same time, only one calculation is needed.
  Vwb.V9 <- data.s3$vol.WB[7] # [m3]
  Awb.V9 <- data.s3$area.WB[7] # [m2]
  veff_static.V9 <- 10^(logKwb$logKwb) * Vwb.V9 * 
    (1 - exp(-ko$ko * Awb.V9 / Vwb.V9 / 10^(logKwb$logKwb) * data.s3$time[7] / 24))
  # 10, 11 and 12 same office, 2 static samples deployed, but different time and size of WB.
  Vwb.1.V10.11.12 <- data.s3$vol.WB[10] # [m3]
  Awb.1.V10.11.12 <- data.s3$area.WB[10] # [m2]
  veff_static.1.V10.11.12 <- 10^(logKwb$logKwb) * Vwb.1.V10.11.12 * 
    (1 - exp(-ko$ko * Awb.1.V10.11.12 / Vwb.1.V10.11.12 / 10^(logKwb$logKwb) * data.s3$time[10] / 24))
  Vwb.2.V10.11.12 <- data.s3$vol.WB[11] # [m3]
  Awb.2.V10.11.12 <- data.s3$area.WB[11] # [m2]
  veff_static.2.V10.11.12 <- 10^(logKwb$logKwb) * Vwb.2.V10.11.12 * 
    (1 - exp(-ko$ko * Awb.2.V10.11.12 / Vwb.2.V10.11.12 / 10^(logKwb$logKwb) * data.s3$time[11] / 24))
  # Volunteers 13 and 14 same office
  Vwb.V13.14 <- data.s3$vol.WB[12] # [m3]
  Awb.V13.14 <- data.s3$area.WB[12] # [m2]
  veff_static.V13.14 <- 10^(logKwb$logKwb) * Vwb.V13.14 * 
    (1 - exp(-ko$ko * Awb.V13.14 / Vwb.V13.14 / 10^(logKwb$logKwb) * data.s3$time[12] / 24))
  
  # Calculate air concentration in ng/m3 from static WBs
  # For V4
  conc.V4 <- as.data.frame(t(data.s3.1[1, ] / veff_static.V4))
  # For V5
  conc.V5 <- as.data.frame(t(data.s3.1[1, ] / veff_static.V5))
  # For V6
  conc.V6 <- as.data.frame(t(data.s3.1[2, ] / veff_static.V6))
  # For V7
  conc.V7.1 <- as.data.frame(t(data.s3.1[3, ] / veff_static.V7.1))
  conc.V7.2 <- as.data.frame(t(data.s3.1[4, ] / veff_static.V7.2))
  conc.V7.3 <- as.data.frame(t(data.s3.1[5, ] / veff_static.V7.3))
  conc.V7 <- (conc.V7.1 + conc.V7.2 + conc.V7.3) / 3
  # For V8
  conc.V8 <- as.data.frame(t(data.s3.1[6, ] / veff_static.V8))
  # For V9 (3 replicates, but the same veff)
  conc.V9 <- data.frame(colMeans(data.s3.1[7:9, ])) / veff_static.V9
  # For V10.11.12 (1)
  conc.1.V10.11.12 <- as.data.frame(t(data.s3.1[10, ] / veff_static.1.V10.11.12))
  # For V10.11.12 (2)
  conc.2.V10.11.12 <- as.data.frame(t(data.s3.1[11, ] / veff_static.2.V10.11.12))
  # For V10, 11 and 12
  conc.V10.11.12 <- (conc.1.V10.11.12 + conc.2.V10.11.12) / 2
  # For V10
  conc.V10 <- conc.V10.11.12
  # For V11
  conc.V11 <- conc.V10.11.12
  # For V12
  conc.V12 <- conc.V10.11.12
  # For V13
  conc.V13 <- as.data.frame(t(data.s3.1[12, ] / veff_static.V13.14))
  # For V14
  conc.V14 <- as.data.frame(t(data.s3.1[12, ] / veff_static.V13.14))
  
  # Combine concentrations
  conc.air <- cbind(conc.V4, conc.V5, conc.V6, conc.V7, conc.V8, conc.V9,
                    conc.V10, conc.V11, conc.V12, conc.V13, conc.V14)
  # Change column names of the last three columns
  colnames(conc.air) <- c("Conc.Air.V4", "Conc.Air.V5", "Conc.Air.V6",
                          "Conc.Air.V7", "Conc.Air.V8", "Conc.Air.V9",
                          "Conc.Air.V10", "Conc.Air.V11", "Conc.Air.V12",
                          "Conc.Air.V13", "Conc.Air.V14")
}

# Check total PCB
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)
# See
tPCB.conc.air

# Export results
write.csv(conc.air,
          file = "Output/Data/Study3/VolunteerConcStaticWBStudy3.csv")

# Estimate air concentration from volunteers WBs --------------------------
# Read ko from PersonalSamplingRates.R
ko.p <- read.csv("Output/Data/Studies1_2/Study2/PersonalAveSR.csv")
ko.p <- ko.p[c(1,7)]

{
  Vwb.V4.nd <- data.s3$vol.WB[13] # [m3]
  Awb.V4.nd <- data.s3$area.WB[13] # [m2]
  veff.V4.nd <- 10^(logKwb$logKwb) * Vwb.V4.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4.nd / Vwb.V4.nd / 10^(logKwb$logKwb) * data.s3$time[13] / 24))
  Vwb.V4.d <- data.s3$vol.WB[14] # [m3]
  Awb.V4.d <- data.s3$area.WB[14] # [m2]
  veff.V4.d <- 10^(logKwb$logKwb) * Vwb.V4.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4.d / Vwb.V4.d / 10^(logKwb$logKwb) * data.s3$time[14] / 24))
  Vwb.V5.nd <- data.s3$vol.WB[15] # [m3]
  Awb.V5.nd <- data.s3$area.WB[15] # [m2]
  veff.V5.nd <- 10^(logKwb$logKwb) * Vwb.V5.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5.nd / Vwb.V5.nd / 10^(logKwb$logKwb) * data.s3$time[15] / 24))
  Vwb.V5.d <- data.s3$vol.WB[16] # [m3]
  Awb.V5.d <- data.s3$area.WB[16] # [m2]
  veff.V5.d <- 10^(logKwb$logKwb) * Vwb.V5.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5.d / Vwb.V5.d / 10^(logKwb$logKwb) * data.s3$time[16] / 24))
  Vwb.V6.nd <- data.s3$vol.WB[17] # [m3]
  Awb.V6.nd <- data.s3$area.WB[17] # [m2]
  veff.V6.nd <- 10^(logKwb$logKwb) * Vwb.V6.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V6.nd / Vwb.V6.nd / 10^(logKwb$logKwb) * data.s3$time[17] / 24))
  Vwb.V6.d <- data.s3$vol.WB[18] # [m3]
  Awb.V6.d <- data.s3$area.WB[18] # [m2]
  veff.V6.d <- 10^(logKwb$logKwb) * Vwb.V6.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V6.d / Vwb.V6.d / 10^(logKwb$logKwb) * data.s3$time[18] / 24))
  Vwb.V7.nd <- data.s3$vol.WB[19] # [m3]
  Awb.V7.nd <- data.s3$area.WB[19] # [m2]
  veff.V7.nd <- 10^(logKwb$logKwb) * Vwb.V7.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V7.nd / Vwb.V7.nd / 10^(logKwb$logKwb) * data.s3$time[19] / 24))
  Vwb.V7.d <- data.s3$vol.WB[20] # [m3]
  Awb.V7.d <- data.s3$area.WB[20] # [m2]
  veff.V7.d <- 10^(logKwb$logKwb) * Vwb.V7.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V7.d / Vwb.V7.d / 10^(logKwb$logKwb) * data.s3$time[20] / 24))
  Vwb.V8.nd <- data.s3$vol.WB[21] # [m3]
  Awb.V8.nd <- data.s3$area.WB[21] # [m2]
  veff.V8.nd <- 10^(logKwb$logKwb) * Vwb.V8.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V8.nd / Vwb.V8.nd / 10^(logKwb$logKwb) * data.s3$time[21] / 24))
  Vwb.V8.d <- data.s3$vol.WB[22] # [m3]
  Awb.V8.d <- data.s3$area.WB[22] # [m2]
  veff.V8.d <- 10^(logKwb$logKwb) * Vwb.V8.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V8.d / Vwb.V8.d / 10^(logKwb$logKwb) * data.s3$time[22] / 24))
  Vwb.V9.nd <- data.s3$vol.WB[23] # [m3]
  Awb.V9.nd <- data.s3$area.WB[23] # [m2]
  veff.V9.nd <- 10^(logKwb$logKwb) * Vwb.V9.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V9.nd / Vwb.V9.nd / 10^(logKwb$logKwb) * data.s3$time[23] / 24))
  Vwb.V9.d <- data.s3$vol.WB[24] # [m3]
  Awb.V9.d <- data.s3$area.WB[24] # [m2]
  veff.V9.d <- 10^(logKwb$logKwb) * Vwb.V9.d * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V9.d / Vwb.V9.d / 10^(logKwb$logKwb) * data.s3$time[24] / 24))
  Vwb.V10.nd <- data.s3$vol.WB[25] # [m3]
  Awb.V10.nd <- data.s3$area.WB[25] # [m2]
  veff.V10.nd <- 10^(logKwb$logKwb) * Vwb.V10.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V10.nd / Vwb.V10.nd / 10^(logKwb$logKwb) * data.s3$time[25] / 24))
  Vwb.V11.nd <- data.s3$vol.WB[26] # [m3]
  Awb.V11.nd <- data.s3$area.WB[26] # [m2]
  veff.V11.nd <- 10^(logKwb$logKwb) * Vwb.V11.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V11.nd / Vwb.V11.nd / 10^(logKwb$logKwb) * data.s3$time[26] / 24))
  Vwb.V12.nd <- data.s3$vol.WB[27] # [m3]
  Awb.V12.nd <- data.s3$area.WB[27] # [m2]
  veff.V12.nd <- 10^(logKwb$logKwb) * Vwb.V12.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V12.nd / Vwb.V12.nd / 10^(logKwb$logKwb) * data.s3$time[27] / 24))
  Vwb.V13.nd <- data.s3$vol.WB[28] # [m3]
  Awb.V13.nd <- data.s3$area.WB[28] # [m2]
  veff.V13.nd <- 10^(logKwb$logKwb) * Vwb.V13.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V13.nd / Vwb.V13.nd / 10^(logKwb$logKwb) * data.s3$time[28] / 24))
  Vwb.V14.nd <- data.s3$vol.WB[29] # [m3]
  Awb.V14.nd <- data.s3$area.WB[29] # [m2]
  veff.V14.nd <- 10^(logKwb$logKwb) * Vwb.V14.nd * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V14.nd / Vwb.V14.nd / 10^(logKwb$logKwb) * data.s3$time[29] / 24))
}

# Estimate air concentration in ng/m3 from WBs
# Using Veff
{
  conc.V4.nd <- as.data.frame(t(data.s3.1[13, ] / veff.V4.nd))
  conc.V4.d <- as.data.frame(t(data.s3.1[14, ] / veff.V4.d))
  conc.V5.nd <- as.data.frame(t(data.s3.1[15, ] / veff.V5.nd))
  conc.V5.d <- as.data.frame(t(data.s3.1[16, ] / veff.V5.d))
  conc.V6.nd <- as.data.frame(t(data.s3.1[17, ] / veff.V6.nd))
  conc.V6.d <- as.data.frame(t(data.s3.1[18, ] / veff.V6.d))
  conc.V7.nd <- as.data.frame(t(data.s3.1[19, ] / veff.V7.nd))
  conc.V7.d <- as.data.frame(t(data.s3.1[20, ] / veff.V7.d))
  conc.V8.nd <- as.data.frame(t(data.s3.1[21, ] / veff.V8.nd))
  conc.V8.d <- as.data.frame(t(data.s3.1[22, ] / veff.V8.d))
  conc.V9.nd <- as.data.frame(t(data.s3.1[23, ] / veff.V9.nd))
  conc.V9.d <- as.data.frame(t(data.s3.1[24, ] / veff.V9.d))
  conc.V10.nd <- as.data.frame(t(data.s3.1[25, ] / veff.V10.nd))
  conc.V11.nd <- as.data.frame(t(data.s3.1[26, ] / veff.V11.nd))
  conc.V12.nd <- as.data.frame(t(data.s3.1[27, ] / veff.V12.nd))
  conc.V13.nd <- as.data.frame(t(data.s3.1[28, ] / veff.V13.nd))
  conc.V14.nd <- as.data.frame(t(data.s3.1[29, ] / veff.V14.nd))
}

# Combined conc WBs
conc.wb <- cbind(conc.V4.nd, conc.V4.d, conc.V5.nd, conc.V5.d, conc.V6.nd, conc.V6.d,
               conc.V7.nd, conc.V7.d, conc.V8.nd, conc.V8.d, conc.V9.nd, conc.V9.d,
               conc.V10.nd, conc.V11.nd, conc.V12.nd, conc.V13.nd, conc.V14.nd)
colnames(conc.wb) <- c("Conc.WB.V4.nd", "Conc.WB.V4.d", "Conc.WB.V5.nd",
                       "Conc.WB.V5.d", "Conc.WB.V6.nd", "Conc.WB.V6.d",
                       "Conc.WB.V7.nd", "Conc.WB.V7.d", "Conc.WB.V8.nd",
                       "Conc.WB.V8.d", "Conc.WB.V9.nd", "Conc.WB.V9.d",
                       "Conc.WB.V10.nd", "Conc.WB.V11.nd", "Conc.WB.V12.nd",
                       "Conc.WB.V13.nd", "Conc.WB.V14.nd")

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- colSums(conc.wb, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Export results
write.csv(conc.wb,
          file = "Output/Data/Study3/VolunteerConcWBStudy3.csv")

# Total PCB plots ---------------------------------------------------------
# Create a data frame with the combined data
data.conc <- data.frame(
  Wb_Concentration = tPCB.conc.wb,
  Air_Concentration = rep(tPCB.conc.air, times = c(2,2,2,2,2,2,1,1,1,1,1)),
  Volunteer = names(tPCB.conc.wb)
)

# Change names for legend
data.conc$Volunteer2 <- data.conc$Volunteer %>%
  str_replace("Conc\\.WB\\.V", "V") %>%
  str_replace("\\.", " ")

# To organize them
data.conc$Volunteer2 <- factor(
  data.conc$Volunteer2,
  levels = str_sort(unique(data.conc$Volunteer2), numeric = TRUE))

color_palette <- c("#377eb8", "#377eb8", "#ff7f0e", "#ff7f0e", "#2ca02c",
                   "#2ca02c", "#d62728", "#d62728", "#9467bd", "#9467bd",
                   "#8c564b", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
                   "#17becf", "black")

shape_palette <- c(21, 24, 21, 24, 21, 24, 21, 24, 21, 24, 21, 24, 24, 24,
                   24, 24, 24)

# Create the plot
plotAirWBtPCB <- ggplot(data.conc, aes(x = Air_Concentration, y = Wb_Concentration,
                                       fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 3.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(1, 10^3), breaks = 10^(0:3),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^3), breaks = 10^(0:3),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  theme(
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right") +
  annotate("text", x = 1, y = 1000, label = "(a)", hjust = 0, vjust = 1,
           size = 6) +
  xlab(expression(bold("Office Air Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted Concentration " *Sigma*"PCB (ng/m"^3*")")))

# Print the plot
print(plotAirWBtPCB)

# Save plot in folder
ggsave("Output/Plots/Study3/VolunteerAirWBtPCBStudy3.png",
       plot = plotAirWBtPCB, width = 6, height = 5, dpi = 500)

# Estimate error (factor of 2) --------------------------------------------
# below 1:1 → model underestimates WB
# above 1:1 → model overestimates WB
# within factor 2 → acceptable agreement
factor_summary <- data.conc %>%
  mutate(factor2 = signif(Wb_Concentration / Air_Concentration, 2)) %>%
  summarise(
    n = n(),
    pct_within_2 = mean(factor2 >= 0.5 & factor2 <= 2, na.rm = TRUE) * 100,
    pct_below_1to1 = mean(factor2 < 1, na.rm = TRUE) * 100,
    pct_above_1to1 = mean(factor2 > 1, na.rm = TRUE) * 100)

factor_summary

# Estimate percentage error ---------------------------------------------------
# Calculate percentage error
percentage_error <-function(observed, predicted) {
  return(abs(observed - predicted)/abs(observed) * 100)
}

# Calculate percentage errors
air_vs_wb_error <- percentage_error(data.conc$Air_Concentration, data.conc$Wb_Concentration)
min(air_vs_wb_error)
max(air_vs_wb_error)
# Calculate mean percent error
mean_error <- mean(air_vs_wb_error)
print(paste("Mean Error:", mean_error))

# Calculate percentage errors btw nd and d WBs (only samples with both nd and d)
# Extract the Wb_Concentration column as a vector
Wb_Concentration <- data.conc[1:12, ]$Wb_Concentration
# Convert the remaining vector into a 2x6 matrix
Wb_matrix <- matrix(Wb_Concentration, nrow = 6, ncol = 2, byrow = TRUE)
colnames(Wb_matrix) <- c('WB_nd', 'WB_d')
# Calculate percentage errors
percentage_error_hand <- percentage_error(Wb_matrix[, 1], Wb_matrix[, 2])
mean(percentage_error_hand)

# Individual PCB Congeners ------------------------------------------------
# Create a data frame with the combined data
# Add PCB as a column
conc_air_long <- conc.air %>%
  rownames_to_column(var = "PCB") %>%
  pivot_longer(
    cols = -PCB,
    names_to = "Volunteer",
    values_to = "Conc.Air"
  ) %>%
  mutate(
    Volunteer = str_extract(Volunteer, "V[0-9]+"))

conc_wb_long <- conc.wb %>%
  rownames_to_column(var = "PCB") %>%
  pivot_longer(
    cols = -PCB,
    names_to = "Volunteer_full",
    values_to = "Conc.WB"
  ) %>%
  mutate(
    Volunteer = str_extract(Volunteer_full, "V[0-9]+"),  # match AIR
    Condition = ifelse(grepl("nd", Volunteer_full), "nd", "d"))

merged_data <- conc_air_long %>%
  inner_join(conc_wb_long, by = c("PCB", "Volunteer"))

filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
# below 1:1 → model underestimates WB
# above 1:1 → model overestimates WB
# within factor 2 → acceptable agreement
factor_summary_congeners <- filtered_data %>%
  mutate(factor2 = signif(Conc.WB / Conc.Air, 2)) %>%
  summarise(
    n = n(),
    pct_within_2 = mean(factor2 >= 0.5 & factor2 <= 2, na.rm = TRUE) * 100,
    pct_below_1to1 = mean(factor2 < 1, na.rm = TRUE) * 100,
    pct_above_1to1 = mean(factor2 > 1, na.rm = TRUE) * 100)
factor_summary_congeners

# Calculate percentage error ---------------------------------------------------
percentage_error <-function(observed, predicted) {
  return(abs(observed - predicted)/abs(observed) * 100)
}

# Calculate percentage errors
percentage_error_pcbi <- percentage_error(filtered_data$Conc.Air, filtered_data$Conc.WB)

# Calculate mean percent error
mean_error <- mean(percentage_error_pcbi)
print(paste("Mean Error:", mean_error))
min_error <- min(percentage_error_pcbi)
print(paste("Minimun Error:", min_error))
max_error <- max(percentage_error_pcbi)
print(paste("Max Error:", max_error))

# Color and shapes from tPCB plot
filtered_data <- filtered_data %>%
  mutate(
    Volunteer2 = paste0("V", gsub("V", "", Volunteer), " ", Condition))

filtered_data$Volunteer2 <- factor(
  filtered_data$Volunteer2,
  levels = c(
    "V4 nd", "V4 d", "V5 nd", "V5 d", "V6 nd", "V6 d",
    "V7 nd", "V7 d", "V8 nd", "V8 d", "V9 nd", "V9 d",
    "V10 nd", "V11 nd", "V12 nd", "V13 nd", "V14 nd"))

color_palette <- setNames(color_palette, levels(filtered_data$Volunteer2))
shape_palette <- setNames(shape_palette, levels(filtered_data$Volunteer2))

# Plot with different colors for each volunteer
plotAirWB <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB,
                                       fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 2.5, color = "black", stroke = 0.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  scale_y_log10(
    limits = c(1e-5, 1e3),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(
    limits = c(1e-5, 1e3),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, color = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, color = "blue", linewidth = 0.7) +
  xlab(expression(bold("Office Air Concentration PCBi (ng/m"^3 * ")"))) +
  ylab(expression(bold("WB Concentration PCBi (ng/m"^3 * ")"))) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 10),
    legend.position = "right") +
  labs(fill = "Volunteers", shape = "Volunteers") +
  annotate("text", x = 10^-5, y = 1000, label = "(b)", hjust = 0, vjust = 1,
           size = 6)

# Print the plot
print(plotAirWB)

# Save plot in folder
ggsave("Output/Plots/Study3/VolunteerAirWBPCBiStudy3.png", plot = plotAirWB, width = 6,
       height = 5, dpi = 500)

# Two plots (facet_wrap)
filtered_data$Condition <- factor(
  filtered_data$Condition,
  levels = c("nd", "d"),
  labels = c("Non-dominant WB", "Dominant WB"))

filtered_data$Volunteer <- factor(filtered_data$Volunteer)

filtered_data$Volunteer <- factor(
  filtered_data$Volunteer,
  levels = paste0("V", 4:14))

vol_levels <- levels(filtered_data$Volunteer)

color_palette <- setNames(
  scales::hue_pal()(length(vol_levels)),
  vol_levels)

shape_palette <- setNames(
  rep(c(21, 22, 23, 24, 25), length.out = length(vol_levels)),
  vol_levels)

plot.2 <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB,
                                       fill = Volunteer, shape = Volunteer,
                                    color = Volunteer)) +
  geom_point(size = 2.5, color = "black", stroke = 0.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  scale_color_manual(values = color_palette) +
  scale_y_log10(
    limits = c(1e-5, 1e3),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(
    limits = c(1e-5, 1e3),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl") +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, color = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, color = "blue", linewidth = 0.7) +
  facet_wrap(~Condition, ncol = 2) +
  xlab(expression(bold("Office Air Concentration PCBi (ng/m"^3 * ")"))) +
  ylab(expression(bold("WB Concentration PCBi (ng/m"^3 * ")"))) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 10),
    legend.position = "right") +
  labs(fill = "Volunteers", shape = "Volunteers")

# Print the plot
print(plot.2)

# Save plot in folder
ggsave("Output/Plots/Study3/VolunteerAirWBPCBi2plotsStudy3.png",
       plot = plot.2, width = 8, height = 5, dpi = 500)
