## Concentration estimation
# Office only

# Install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
}

# Read measured values from excel -----------------------------------------
{
  data <- data.frame(read_excel("Data/Volunteers.xlsx", sheet = "Sheet1",
                                col_names = TRUE, col_types = NULL))
  data.2 <- data.frame(read_excel("Data/VolunteersV02.xlsx", sheet = "Sheet1",
                                  col_names = TRUE, col_types = NULL))
  logKoa <- data.frame(read_excel("Data/logKoa.xlsx", sheet = "logKoa",
                                  col_names = TRUE, col_types = NULL))
  # ko from SamplingRates_ko.R file
  ko <- read.csv("Output/Data/csv/SamplingRates/SR/WDSamplingRateStatV3.csv")
  # Select only ko [m/d]
  ko <- ko[c(2,6)]
}

# Calculate logKws
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
logKwb <- data.frame(
  congener = logKoa$congener,
  logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96

# Volunteers 1 to 8 (V)
# Calculate air PCB concentration from static WBs -------------------------
{
  # Calculate effective volume for static WBs
  # Use effective volume. Adult WBs
  Vwb <- 4.73 * 10^-6 # [m3]
  Awb <- 0.0054773 # [m2]
  veff_static.V1.V2 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[8]))
  veff_static.V3 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[7]))
  veff_static.V4 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[14]))
  veff_static.V5 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[13]))
  veff_static.V6 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[19]))
  veff_static.V7 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[22]))
  veff_static.V8 <- 10^(logKwb$logKwb) * Vwb * 
    (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.2$office.time.day[3]))
  
  # Calculate air concentration in ng/m3 from static WBs
  # For V1 & V2
  conc.V1.V2 <- as.data.frame(t(data[8, 4:176] / veff_static.V1.V2))
  # For V3
  conc.V3 <- as.data.frame(t(data[7, 4:176] / veff_static.V3))
  # For V4
  conc.V4 <- data.frame(colMeans(data[14:15, 4:176])) / veff_static.V4
  # For V5
  conc.V5 <- as.data.frame(t(data[13, 4:176] / veff_static.V5))
  # For V6
  conc.V6 <- data.frame(colMeans(data[19:21, 4:176])) / veff_static.V6
  # For V7
  conc.V7 <- data.frame(colMeans(data[22:23, 4:176])) / veff_static.V7
  # For V8
  conc.V8 <- as.data.frame(t(data.2[3, 4:176] / veff_static.V8))
  
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


# Estimate air concentration from volunteers WBs --------------------------
Vwb <- 4.73 * 10^-6 # [m3]
Awb <- 0.0054773 # [m2]
veff.V1.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[1]))
veff.V1.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[2]))
veff.V2.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[3]))
veff.V2.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[4]))
veff.V3.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[5]))
veff.V3.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[6]))
veff.V4.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[9]))
veff.V4.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[10]))
veff.V5.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[11]))
veff.V5.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[12]))
veff.V6.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[17]))
veff.V6.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[18]))
veff.V7.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[25]))
veff.V7.r <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data$time.day[24]))
veff.V8.l <- 10^(logKwb$logKwb) * Vwb * 
  (1 - exp(-ko$ko * Awb / Vwb / 10^(logKwb$logKwb) * data.2$office.time.day[13]))

# Estimate air concentration in ng/m3 from WBs
# For V1
conc.V1.l <- as.data.frame(t(data[1, 4:176] / veff.V1.l))
conc.V1.r <- as.data.frame(t(data[2, 4:176] / veff.V1.r))
conc.V2.l <- as.data.frame(t(data[3, 4:176] / veff.V2.l))
conc.V2.r <- as.data.frame(t(data[4, 4:176] / veff.V2.l))
conc.V3.l <- as.data.frame(t(data[5, 4:176] / veff.V3.l))
conc.V3.r <- as.data.frame(t(data[6, 4:176] / veff.V3.r))
conc.V4.l <- as.data.frame(t(data[9, 4:176] / veff.V4.l))
conc.V4.r <- as.data.frame(t(data[10, 4:176] / veff.V4.r))
conc.V5.l <- as.data.frame(t(data[11, 4:176] / veff.V5.l))
conc.V5.r <- as.data.frame(t(data[12, 4:176] / veff.V5.r))
conc.V6.l <- as.data.frame(t(data[17, 4:176] / veff.V6.l))
conc.V6.r <- as.data.frame(t(data[18, 4:176] / veff.V6.r))
conc.V7.l <- as.data.frame(t(data[25, 4:176] / veff.V7.l))
conc.V7.r <- as.data.frame(t(data[24, 4:176] / veff.V7.r))
conc.V8.l <- as.data.frame(t(data.2[13, 4:176] / veff.V8.l))





# Combined wore WBs
wb.wr <- rbind(conc.V1.l, conc.V1.r, wb.Ya.l, wb.Ya.r, wb.Ea.l, wb.Ea.r,
               wb.Cr.l, wb.Cr.r, wb.Hu.l, wb.Hu.r, wb.Xu.l, wb.Xu.r,
               wb.Gi.l, wb.Gi.r, wb.Xue.l)

# Estimate air concentration using SR, mass of wore WBs & time ------------
# Extract congener names from sr
congener_names <- sr$congener
# Subset wb.wr to include only the common congeners
wb_common <- wb.wr[, intersect(colnames(wb.wr), congener_names)]
# Subset sr to include only the common congeners
sr_common <- sr[sr$congener %in% colnames(wb_common), ]
# Ensure that sr_common is ordered to match the columns of wb_common
sr_common <- sr_common[match(colnames(wb_common), sr_common$congener), ]
# Divide each element in wb_common by corresponding element in sr_common$Average_Sampling_Rate
wb_div_sr <- sweep(wb_common, 2, sr_common$Average_Sampling_Rate, FUN = "/")
# Alternative method -> use a constant sampling rate
# wb_div_sr <- sweep(wb_common, 2, 1.5, FUN = "/")
# Extract time.day from wb.wr
wb_time_day <- wb.wr$time.day
# Divide wb_div_sr further by corresponding value in wb_time_day
conc.wb <- sweep(wb_div_sr, 1, wb_time_day, FUN = "/")
rownames(conc.wb) <- c('wb.Mi.l', 'wb.Mi.r', 'wb.Ya.l', 'wb.Ya.r', 'wb.Ea.l',
                       'wb.Ea.r', 'wb.Cr.l', 'wb.Cr.r', 'wb.Hu.l', 'wb.Hu.r',
                       'wb.Xu.l', 'wb.Xu.r', 'wb.Gi.l', 'wb.Gi.r', 'wb.Xue.l')

# Match congeners in both dataset -----------------------------------------
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb), rownames(conc.air))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.air))
# Subset conc.air to include only the rows with matching row names
conc_air_common <- conc.air[matching_indices, ]

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- rowSums(conc.wb, na.rm = TRUE)
# Convert list column to numeric, replacing NA values with 0
conc_air_common_numeric <- as.data.frame(sapply(conc_air_common, as.numeric))
tPCB.conc.air <- colSums(conc_air_common_numeric, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Plots -------------------------------------------------------------------
# tPCB
# Create a data frame with the combined data
data <- data.frame(
  Wb_Concentration = tPCB.conc.wb,
  Air_Concentration = rep(tPCB.conc.air, times = c(4, 2, 2, 2, 2, 2, 1)),
  Volunteer = names(tPCB.conc.wb)
)

# Change names for legend
data$Volunteer2 <- case_when(
  data$Volunteer == "wb.Mi.l" ~ "Vol. 1 nd",
  data$Volunteer == "wb.Mi.r" ~ "Vol. 1 d",
  data$Volunteer == "wb.Ya.l" ~ "Vol. 2 d",
  data$Volunteer == "wb.Ya.r" ~ "Vol. 2 nd",
  data$Volunteer == "wb.Ea.l" ~ "Vol. 3 nd",
  data$Volunteer == "wb.Ea.r" ~ "Vol. 3 d",
  data$Volunteer == "wb.Cr.l" ~ "Vol. 4 nd",
  data$Volunteer == "wb.Cr.r" ~ "Vol. 4 d",
  data$Volunteer == "wb.Hu.l" ~ "Vol. 5 d",
  data$Volunteer == "wb.Hu.r" ~ "Vol. 5 nd",
  data$Volunteer == "wb.Xu.l" ~ "Vol. 6 nd",
  data$Volunteer == "wb.Xu.r" ~ "Vol. 6 d",
  data$Volunteer == "wb.Gi.l" ~ "Vol. 7 nd",
  data$Volunteer == "wb.Gi.r" ~ "Vol. 7 d",
  data$Volunteer == "wb.Xue.l" ~ "Vol. 8 nd",
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
plotAirWBtPCB <- ggplot(data, aes(x = Air_Concentration, y = Wb_Concentration,
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
        axis.title.x = element_text(face = "bold", size = 14)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  theme(legend.position = "right") # Optional: Adjust legend position

# Print the plot
print(plotAirWBtPCB)

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VolunteerAirWBtPCBV2.png", plot = plotAirWBtPCB, width = 6,
       height = 5, dpi = 500)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
data$factor2 <- data$Wb_Concentration/data$Air_Concentration

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(data[data$factor2 > 0.5 & data$factor2 < 2, ])/nrow(data)*100

# Estimate percentage error ---------------------------------------------------
# Calculate percentage error
percentage_error <-function(observed, predicted) {
  return(abs(observed - predicted)/abs(observed) * 100)
}

# Calculate percentage errors
percentage_error <- percentage_error(data$Air_Concentration, data$Wb_Concentration)
min(percentage_error)
max(percentage_error)
# Calculate mean percent error
mean_error <- mean(percentage_error)
print(paste("Mean Error:", mean_error))

# Calculate percentage errors btw left and right WBs
# Extract the Wb_Concentration column as a vector
Wb_Concentration <- data$Wb_Concentration
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
# Reshape conc_air_common to long format

# Add PCB as a column
conc_air_common <- conc_air_common %>% rownames_to_column(var = "PCB")

conc_air_long <- conc_air_common %>%
  pivot_longer(
    cols = -PCB,  # All columns except 'PCB'
    names_to = "Volunteer",
    values_to = "Conc.Air"
  ) %>%
  mutate(Volunteer = gsub("Conc.Air.", "", Volunteer))

# Duplicate Mi.Ya rows and create new rows for Mi and Ya
conc_air_long <- conc_air_long %>%
  # Filter out Mi.Ya rows
  filter(Volunteer != "Mi.Ya") %>%
  # Create new rows for Mi and Ya
  bind_rows(
    conc_air_long %>%
      filter(Volunteer == "Mi.Ya") %>%
      # Duplicate the rows and assign Mi and Ya
      uncount(2) %>%
      mutate(Volunteer = ifelse(row_number() %% 2 == 1, "Mi", "Ya"))
  )

conc_wb_long <- as.data.frame(conc.wb) %>%
  rownames_to_column(var = "Volunteer") %>%  # Convert row names to a column
  pivot_longer(
    cols = -Volunteer,  # All columns except 'Volunteer'
    names_to = "PCB",   # Column for PCB names
    values_to = "Conc.WB"  # Column for concentration values
  ) %>%
  mutate(Volunteer = gsub("wb.", "", Volunteer)) %>%  # Remove 'wb.' from Volunteer names
  select(PCB, Conc.WB, Volunteer)

# Create a simplified Volunteer column in conc_wb_long
conc_wb_long <- conc_wb_long %>%
  mutate(Volunteer_Simplified = gsub("\\..*", "", Volunteer))  # Remove everything after and including "."

# Merge the two data frames
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
percentage_error <- percentage_error(filtered_data$Conc.Air, filtered_data$Conc.WB)

# Calculate mean percent error
mean_error <- mean(percentage_error)
print(paste("Mean Error:", mean_error))
min_error <- min(percentage_error)
print(paste("Minimun Error:", min_error))
max_error <- max(percentage_error)
print(paste("Max Error:", max_error))

# Define the threshold
threshold <- 200

# Count how many values are above the threshold
count_above_threshold <- sum(percentage_error > threshold)

# Get the total number of values
total_values <- length(percentage_error)

# Calculate the proportion
proportion_above_threshold <- count_above_threshold / total_values

# Print the results
cat("Number of values above", threshold, ":", count_above_threshold, "\n")
cat("Total number of values:", total_values, "\n")
cat("Proportion of values above", threshold, ":", proportion_above_threshold, "\n")

# Color and shapes from tPCB plot
# Change names for legend
filtered_data$Volunteer.y <- case_when(
  filtered_data$Volunteer.y == "Mi.l" ~ "Vol. 1 nd",
  filtered_data$Volunteer.y == "Mi.r" ~ "Vol. 1 d",
  filtered_data$Volunteer.y == "Ya.l" ~ "Vol. 2 d",
  filtered_data$Volunteer.y == "Ya.r" ~ "Vol. 2 nd",
  filtered_data$Volunteer.y == "Ea.l" ~ "Vol. 3 nd",
  filtered_data$Volunteer.y == "Ea.r" ~ "Vol. 3 d",
  filtered_data$Volunteer.y == "Cr.l" ~ "Vol. 4 nd",
  filtered_data$Volunteer.y == "Cr.r" ~ "Vol. 4 d",
  filtered_data$Volunteer.y == "Hu.l" ~ "Vol. 5 d",
  filtered_data$Volunteer.y == "Hu.r" ~ "Vol. 5 nd",
  filtered_data$Volunteer.y == "Xu.l" ~ "Vol. 6 nd",
  filtered_data$Volunteer.y == "Xu.r" ~ "Vol. 6 d",
  filtered_data$Volunteer.y == "Gi.l" ~ "Vol. 7 nd",
  filtered_data$Volunteer.y == "Gi.r" ~ "Vol. 7 d",
  filtered_data$Volunteer.y == "Xue.l" ~ "Vol. 8 nd",
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
          axis.title.x = element_text(face = "bold", size = 14)) +
    geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
    geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
    geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
    guides(fill = guide_legend(override.aes = list(color = NA))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  theme(legend.position = "right")

# Print the plots
print(plotAirWBPCBi)

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VolunteerAirWBPCBiV2.png", plot = plotAirWBPCBi, width = 6,
       height = 5, dpi = 500)


