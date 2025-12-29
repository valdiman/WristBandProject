## Concentration estimation
# Office only

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
}

# Read data ---------------------------------------------------------------
{
  data.0 <- read.csv("Data/IRO/SampleMassStudy3_4_5.csv")
  data <- data.0[1:25, c(1, 7, 9:10, 12:184)]
  data.2 <- data.0[26:38, c(1, 7, 9:10, 12:184)]
  logKoa <- read.csv("Data/logKoa.csv")
  # ko from SamplingRates_ko.R file
  ko <- read.csv("Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV1.csv")
  # Select only ko [m/d]
  ko <- ko[c(2,6)]
}

# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# Select only congeners
common_ids <- intersect(ko$congener, logKwb$congener)
ko <- ko[ko$congener %in% common_ids, ]
logKwb <- logKwb[logKwb$congener %in% common_ids, ]
data.1 <- data[, common_ids]
data.2.1 <- data.2[, common_ids]

# Volunteers 1 to 8 (V)
# Calculate air PCB concentration from static WBs -------------------------
{
  # Calculate effective volume for static WBs
  # Use effective volume. Adult WBs
  Vwb.V1.V2 <- data$vol.WB[8] # [m3]
  Awb.V1.V2 <- data$area.WB[8] # [m2]
  veff_static.V1.V2 <- 10^(logKwb$logKwb) * Vwb.V1.V2 * 
    (1 - exp(-ko$ko * Awb.V1.V2 / Vwb.V1.V2 / 10^(logKwb$logKwb) * data$time[8] / 24))
  Vwb.V3 <- data$vol.WB[7] # [m3]
  Awb.V3 <- data$area.WB[7] # [m2]
  veff_static.V3 <- 10^(logKwb$logKwb) * Vwb.V3 * 
    (1 - exp(-ko$ko * Awb.V3 / Vwb.V3 / 10^(logKwb$logKwb) * data$time[7] / 24))
  Vwb.V4 <- data$vol.WB[14] # [m3]
  Awb.V4 <- data$area.WB[14] # [m2]
  veff_static.V4 <- 10^(logKwb$logKwb) * Vwb.V4 * 
    (1 - exp(-ko$ko * Awb.V4 / Vwb.V4 / 10^(logKwb$logKwb) * data$time[14] / 24))
  Vwb.V5 <- data$vol.WB[13] # [m3]
  Awb.V5 <- data$area.WB[13] # [m2]
  veff_static.V5 <- 10^(logKwb$logKwb) * Vwb.V5 * 
    (1 - exp(-ko$ko * Awb.V5 / Vwb.V5 / 10^(logKwb$logKwb) * data$time[13] /24))
  Vwb.V6 <- data$vol.WB[19] # [m3]
  Awb.V6 <- data$area.WB[19] # [m2]
  veff_static.V6 <- 10^(logKwb$logKwb) * Vwb.V6 * 
    (1 - exp(-ko$ko * Awb.V6 / Vwb.V6 / 10^(logKwb$logKwb) * data$time[19] / 24))
  Vwb.V7 <- data$vol.WB[22] # [m3]
  Awb.V7 <- data$area.WB[22] # [m2]
  veff_static.V7 <- 10^(logKwb$logKwb) * Vwb.V7 * 
    (1 - exp(-ko$ko * Awb.V7 / Vwb.V7 / 10^(logKwb$logKwb) * data$time[22] / 24))
  Vwb.V8 <- data.2$vol.WB[3] # [m3]
  Awb.V8 <- data.2$area.WB[3] # [m2]
  veff_static.V8 <- 10^(logKwb$logKwb) * Vwb.V8 * 
    (1 - exp(-ko$ko * Awb.V8 / Vwb.V8 / 10^(logKwb$logKwb) * data.2$time[3] / 24))
  
  # Calculate air concentration in ng/m3 from static WBs
  # For V1 & V2
  conc.V1.V2 <- as.data.frame(t(data.1[8, ] / veff_static.V1.V2))
  # For V3
  conc.V3 <- as.data.frame(t(data.1[7, ] / veff_static.V3))
  # For V4
  conc.V4 <- data.frame(colMeans(data.1[14:16, ])) / veff_static.V4
  # For V5
  conc.V5 <- as.data.frame(t(data.1[13, ] / veff_static.V5))
  # For V6
  conc.V6 <- data.frame(colMeans(data.1[19:21, ])) / veff_static.V6
  # For V7
  conc.V7 <- data.frame(colMeans(data.1[22:23, ])) / veff_static.V7
  # For V8
  conc.V8 <- as.data.frame(t(data.2.1[3, ] / veff_static.V8))
  
  # Combine concentrations
  conc.air <- cbind(conc.V1.V2, conc.V3, conc.V4, conc.V5, conc.V6, conc.V7,
                    conc.V8)
  # Change column names of the last three columns
  colnames(conc.air) <- c("Conc.Air.V1.V2", "Conc.Air.V3", "Conc.Air.V4",
                          "Conc.Air.V5", "Conc.Air.V6", "Conc.Air.V7",
                          "Conc.Air.V8")
}

# Check total PCB
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)
# See
tPCB.conc.air

# Export results
write.csv(conc.air,
          file = "Output/Data/csv/Volunteer/VolunteerConcStaticWB.csv")

# Estimate air concentration from volunteers WBs --------------------------
# Read ko from PersonalSamplingRates
ko.p <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSR.csv")
ko.p <- ko.p[c(1,7)]

{
  Vwb.V1.l <- data$vol.WB[1] # [m3]
  Awb.V1.l <- data$area.WB[1] # [m2]
  veff.V1.l <- 10^(logKwb$logKwb) * Vwb.V1.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V1.l / Vwb.V1.l / 10^(logKwb$logKwb) * data$time[1] / 24))
  Vwb.V1.r <- data$vol.WB[2] # [m3]
  Awb.V1.r <- data$area.WB[2] # [m2]
  veff.V1.r <- 10^(logKwb$logKwb) * Vwb.V1.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V1.r / Vwb.V1.r / 10^(logKwb$logKwb) * data$time[2] / 24))
  Vwb.V2.l <- data$vol.WB[3] # [m3]
  Awb.V2.l <- data$area.WB[3] # [m2]
  veff.V2.l <- 10^(logKwb$logKwb) * Vwb.V2.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V2.l / Vwb.V2.l / 10^(logKwb$logKwb) * data$time[3] / 24))
  Vwb.V2.r <- data$vol.WB[4] # [m3]
  Awb.V2.r <- data$area.WB[4] # [m2]
  veff.V2.r <- 10^(logKwb$logKwb) * Vwb.V2.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V2.r / Vwb.V2.r / 10^(logKwb$logKwb) * data$time[4] / 24))
  Vwb.V3.l <- data$vol.WB[5] # [m3]
  Awb.V3.l <- data$area.WB[5] # [m2]
  veff.V3.l <- 10^(logKwb$logKwb) * Vwb.V3.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V3.l / Vwb.V3.l / 10^(logKwb$logKwb) * data$time[5] / 24))
  Vwb.V3.r <- data$vol.WB[6] # [m3]
  Awb.V3.r <- data$area.WB[6] # [m2]
  veff.V3.r <- 10^(logKwb$logKwb) * Vwb.V3.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V3.r / Vwb.V3.r / 10^(logKwb$logKwb) * data$time[6] / 24))
  Vwb.V4.l <- data$vol.WB[9] # [m3]
  Awb.V4.l <- data$area.WB[9] # [m2]
  veff.V4.l <- 10^(logKwb$logKwb) * Vwb.V4.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4.l / Vwb.V4.l / 10^(logKwb$logKwb) * data$time[9] / 24))
  Vwb.V4.r <- data$vol.WB[10] # [m3]
  Awb.V4.r <- data$area.WB[10] # [m2]
  veff.V4.r <- 10^(logKwb$logKwb) * Vwb.V4.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4.r / Vwb.V4.r / 10^(logKwb$logKwb) * data$time[10] / 24))
  Vwb.V5.l <- data$vol.WB[11] # [m3]
  Awb.V5.l <- data$area.WB[11] # [m2]
  veff.V5.l <- 10^(logKwb$logKwb) * Vwb.V5.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5.l / Vwb.V5.l / 10^(logKwb$logKwb) * data$time[11] / 24))
  Vwb.V5.r <- data$vol.WB[12] # [m3]
  Awb.V5.r <- data$area.WB[12] # [m2]
  veff.V5.r <- 10^(logKwb$logKwb) * Vwb.V5.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5.r / Vwb.V5.r / 10^(logKwb$logKwb) * data$time[12] / 24))
  Vwb.V6.l <- data$vol.WB[17] # [m3]
  Awb.V6.l <- data$area.WB[17] # [m2]
  veff.V6.l <- 10^(logKwb$logKwb) * Vwb.V6.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V6.l / Vwb.V6.l / 10^(logKwb$logKwb) * data$time[17] / 24))
  Vwb.V6.r <- data$vol.WB[18] # [m3]
  Awb.V6.r <- data$area.WB[18] # [m2]
  veff.V6.r <- 10^(logKwb$logKwb) * Vwb.V6.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V6.r / Vwb.V6.r / 10^(logKwb$logKwb) * data$time[18] / 24))
  Vwb.V7.l <- data$vol.WB[25] # [m3]
  Awb.V7.l <- data$area.WB[25] # [m2]
  veff.V7.l <- 10^(logKwb$logKwb) * Vwb.V7.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V7.l / Vwb.V7.l / 10^(logKwb$logKwb) * data$time[25] / 24))
  Vwb.V7.r <- data$vol.WB[24] # [m3]
  Awb.V7.r <- data$area.WB[24] # [m2]
  veff.V7.r <- 10^(logKwb$logKwb) * Vwb.V7.r * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V7.r / Vwb.V7.r / 10^(logKwb$logKwb) * data$time[24] / 24))
  Vwb.V8.l <- data.2$vol.WB[13] # [m3]
  Awb.V8.l <- data.2$area.WB[13] # [m2]
  veff.V8.l <- 10^(logKwb$logKwb) * Vwb.V8.l * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V8.l / Vwb.V8.l / 10^(logKwb$logKwb) * data.2$time[13] / 24))
}

# Estimate air concentration in ng/m3 from WBs
# Using Veff
{
  conc.V1.l <- as.data.frame(t(data.1[1, ] / veff.V1.l))
  conc.V1.r <- as.data.frame(t(data.1[2, ] / veff.V1.r))
  conc.V2.l <- as.data.frame(t(data.1[3, ] / veff.V2.l))
  conc.V2.r <- as.data.frame(t(data.1[4, ] / veff.V2.l))
  conc.V3.l <- as.data.frame(t(data.1[5, ] / veff.V3.l))
  conc.V3.r <- as.data.frame(t(data.1[6, ] / veff.V3.r))
  conc.V4.l <- as.data.frame(t(data.1[9, ] / veff.V4.l))
  conc.V4.r <- as.data.frame(t(data.1[10, ] / veff.V4.r))
  conc.V5.l <- as.data.frame(t(data.1[11, ] / veff.V5.l))
  conc.V5.r <- as.data.frame(t(data.1[12, ] / veff.V5.r))
  conc.V6.l <- as.data.frame(t(data.1[17, ] / veff.V6.l))
  conc.V6.r <- as.data.frame(t(data.1[18, ] / veff.V6.r))
  conc.V7.l <- as.data.frame(t(data.1[25, ] / veff.V7.l))
  conc.V7.r <- as.data.frame(t(data.1[24, ] / veff.V7.r))
  conc.V8.l <- as.data.frame(t(data.2.1[13, ] / veff.V8.l))
}

# Combined conc WBs
conc.wb <- cbind(conc.V1.l, conc.V1.r, conc.V2.l, conc.V2.r, conc.V3.l, conc.V3.r,
               conc.V4.l, conc.V4.r, conc.V5.l, conc.V5.r, conc.V6.l, conc.V6.r,
               conc.V7.l, conc.V7.r, conc.V8.l)
colnames(conc.wb) <- c('conc.V1.l', 'conc.V1.r', 'conc.V2.l', 'conc.V2.r', 'conc.V3.l',
                       'conc.V3.r', 'conc.V4.l', 'conc.V4.r', 'conc.V5.l', 'conc.V5.r',
                       'conc.V6.l', 'conc.V6.r', 'conc.V7.l', 'conc.V7.r', 'conc.V8.l')

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- colSums(conc.wb, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Export results
write.csv(conc.wb,
          file = "Output/Data/csv/Volunteer/VolunteerConcWB.csv")

# Total PCB plots ---------------------------------------------------------
# Create a data frame with the combined data
data.conc <- data.frame(
  Wb_Concentration = tPCB.conc.wb,
  Air_Concentration = rep(tPCB.conc.air, times = c(4, 2, 2, 2, 2, 2, 1)),
  Volunteer = names(tPCB.conc.wb)
)

# Change names for legend
data.conc$Volunteer2 <- case_when(
  data.conc$Volunteer == "conc.V1.l" ~ "Vol. 1 nd",
  data.conc$Volunteer == "conc.V1.r" ~ "Vol. 1 d",
  data.conc$Volunteer == "conc.V2.l" ~ "Vol. 2 d",
  data.conc$Volunteer == "conc.V2.r" ~ "Vol. 2 nd",
  data.conc$Volunteer == "conc.V3.l" ~ "Vol. 3 nd",
  data.conc$Volunteer == "conc.V3.r" ~ "Vol. 3 d",
  data.conc$Volunteer == "conc.V4.l" ~ "Vol. 4 nd",
  data.conc$Volunteer == "conc.V4.r" ~ "Vol. 4 d",
  data.conc$Volunteer == "conc.V5.l" ~ "Vol. 5 d",
  data.conc$Volunteer == "conc.V5.r" ~ "Vol. 5 nd",
  data.conc$Volunteer == "conc.V6.l" ~ "Vol. 6 nd",
  data.conc$Volunteer == "conc.V6.r" ~ "Vol. 6 d",
  data.conc$Volunteer == "conc.V7.l" ~ "Vol. 7 nd",
  data.conc$Volunteer == "conc.V7.r" ~ "Vol. 7 d",
  data.conc$Volunteer == "conc.V8.l" ~ "Vol. 8 nd",
  TRUE ~ NA_character_  # This handles any unmatched cases
)

# Define a color palette with enough distinct colors for the number of volunteers
color_palette <- c(
  "#377eb8", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#e41a1c", "#264c73",
  "#f781bf", "#a65628", "#4b0082")

# Define a shape palette with enough distinct shapes for the number of volunteers
shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22, 25)

# Create the plot
plotAirWBtPCB <- ggplot(data.conc, aes(x = Air_Concentration, y = Wb_Concentration,
                                  fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 3.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(1, 10^3),
                breaks = 10^(0:3),  # Integer powers of 10 only
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^3),
                breaks = 10^(0:3),  # Integer powers of 10 only
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Air Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  theme(legend.position = "right") + 
  annotate("text", x = 1, y = 1000,
           label = "(a)", hjust = 0, vjust = 1, 
           size = 6, color = "black")


# Print the plot
print(plotAirWBtPCB)

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VolunteerAirWBtPCB.png",
       plot = plotAirWBtPCB, width = 6, height = 5, dpi = 500)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
data.conc$factor2 <- data.conc$Wb_Concentration/data.conc$Air_Concentration

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(data.conc[data.conc$factor2 > 0.5 & data.conc$factor2 < 2, ])/nrow(data.conc)*100

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

# Calculate percentage errors btw left and right WBs
# Extract the Wb_Concentration column as a vector
Wb_Concentration <- data.conc$Wb_Concentration
# Remove the last row to ensure an even number of elements
Wb_Concentration <- Wb_Concentration[-length(Wb_Concentration)]
# Convert the remaining vector into a 2x7 matrix
Wb_matrix <- matrix(Wb_Concentration, nrow = 7, ncol = 2, byrow = TRUE)
colnames(Wb_matrix) <- c('WBleft', 'WBright')
# Calculate percentage errors
percentage_error_hand <- percentage_error(Wb_matrix[, 1], Wb_matrix[, 2])
mean(percentage_error_hand)

# Individual PCB Congeners ------------------------------------------------
# Create a data frame with the combined data
# Add PCB as a column
conc_air_long <- conc.air %>% rownames_to_column(var = "PCB")

conc_air_long <- conc_air_long %>%
  pivot_longer(
    cols = -PCB,  # All columns except 'PCB'
    names_to = "Volunteer",
    values_to = "Conc.Air"
  ) %>%
  mutate(Volunteer = gsub("Conc.Air.", "", Volunteer))

# Duplicate Mi.Ya rows and create new rows for Mi and Ya
conc_air_long <- conc_air_long %>%
  # Filter out V1.V2 rows
  filter(Volunteer != "V1.V2") %>%
  # Create new rows for Mi and Ya
  bind_rows(
    conc_air_long %>%
      filter(Volunteer == "V1.V2") %>%
      # Duplicate the rows and assign Mi and Ya
      uncount(2) %>%
      mutate(Volunteer = ifelse(row_number() %% 2 == 1, "V1", "V2"))
  )

conc_wb_long <- as.data.frame(conc.wb) %>%
  rownames_to_column(var = "PCB") %>%  # Convert row names to a column
  pivot_longer(
    cols = -PCB,  # All columns except 'Volunteer'
    names_to = "Volunteer",   # Column for PCB names
    values_to = "Conc.WB"  # Column for concentration values
  ) %>%
  mutate(Volunteer = gsub("wb.", "", Volunteer)) %>%  # Remove 'wb.' from Volunteer names
  select(PCB, Conc.WB, Volunteer)

# Modify Volunteer to remove the 'conc.' prefix and the suffix (.l or .r)
conc_wb_long <- conc_wb_long %>%
  mutate(Volunteer_Simplified = gsub("^conc\\.", "", Volunteer),  # Remove 'conc.' prefix
         Volunteer_Simplified = gsub("\\..*", "", Volunteer_Simplified))  # Remove everything after the first '.'

merged_data <- conc_air_long %>%
  inner_join(conc_wb_long, by = c("PCB", "Volunteer" = "Volunteer_Simplified"))

# Filter out rows with 0 in either Conc.Air or Conc.WB
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
filtered_data$factor2 <- filtered_data$Conc.WB/filtered_data$Conc.Air

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(filtered_data[filtered_data$factor2 > 0.5 & filtered_data$factor2 < 2, ])/nrow(filtered_data)*100

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
# Change names for legend
filtered_data$Volunteer.y <- case_when(
  filtered_data$Volunteer.y == "conc.V1.l" ~ "Vol. 1 nd",
  filtered_data$Volunteer.y == "conc.V1.r" ~ "Vol. 1 d",
  filtered_data$Volunteer.y == "conc.V2.l" ~ "Vol. 2 d",
  filtered_data$Volunteer.y == "conc.V2.r" ~ "Vol. 2 nd",
  filtered_data$Volunteer.y == "conc.V3.l" ~ "Vol. 3 nd",
  filtered_data$Volunteer.y == "conc.V3.r" ~ "Vol. 3 d",
  filtered_data$Volunteer.y == "conc.V4.l" ~ "Vol. 4 nd",
  filtered_data$Volunteer.y == "conc.V4.r" ~ "Vol. 4 d",
  filtered_data$Volunteer.y == "conc.V5.l" ~ "Vol. 5 d",
  filtered_data$Volunteer.y == "conc.V5.r" ~ "Vol. 5 nd",
  filtered_data$Volunteer.y == "conc.V6.l" ~ "Vol. 6 nd",
  filtered_data$Volunteer.y == "conc.V6.r" ~ "Vol. 6 d",
  filtered_data$Volunteer.y == "conc.V7.l" ~ "Vol. 7 nd",
  filtered_data$Volunteer.y == "conc.V7.r" ~ "Vol. 7 d",
  filtered_data$Volunteer.y == "conc.V8.l" ~ "Vol. 8 nd",
  TRUE ~ NA_character_  # Handles any unmatched cases
)

# Plot with different colors for each volunteer  
plotAirWBPCBi <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
                             fill = Volunteer.y, shape = Volunteer.y)) +
    geom_point(size = 2.5, color = "black", stroke = 0.5) +
    theme_bw() +
    theme(aspect.ratio = 15/15) +
    annotation_logticks(sides = "bl") +
    scale_y_log10(limits = c(0.00001, 10^3),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(limits = c(0.00001, 10^3),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(expression(bold("Air Concentration PCBi (ng/m"^3*")"))) +
    ylab(expression(bold("Predicted Concentration PCBi (ng/m"^3*")"))) +
    theme(axis.text.y = element_text(face = "bold", size = 14),
          axis.title.y = element_text(face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 12)) +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
    guides(fill = guide_legend(override.aes = list(color = NA))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  theme(legend.position = "right") +
  annotate("text", x = 0.00001, y = 1000,
           label = "(b)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(plotAirWBPCBi)

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VolunteerAirWBPCBi.png", plot = plotAirWBPCBi, width = 6,
       height = 5, dpi = 500)

