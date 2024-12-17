# Concentration estimation

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
data <- data.frame(read_excel("Data/VolunteersV02.xlsx", sheet = "Sheet1",
                                       col_names = TRUE, col_types = NULL))

# Calculate air PCB concentration from static WBs -------------------------
{
  # Stat samples
  data.Stat.1 <- data.frame(data[1, 3:175])
  data.Stat.2 <- data.frame(data[2, 3:175])
  data.Stat.3 <- data.frame(data[3, 3:175])
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc.air.1 <- as.data.frame(t(data.Stat.1/(0.5*data$time.day[1])))
  conc.air.2 <- as.data.frame(t(data.Stat.1/(0.5*data$time.day[2])))
  conc.air.3 <- as.data.frame(t(data.Stat.1/(0.5*data$time.day[3])))
  # Combine the three conc.air data frames into a list
  conc.air <- list(conc.air.1, conc.air.2, conc.air.3)
  # Bind the data frames together row-wise
  conc.air <- do.call(cbind, conc.air)
  # Calculate the row means for each column and create a new data frame
  conc.air.avg <- as.data.frame(rowMeans(conc.air))
  
  # Change column names of the last three columns
  colnames(conc.air.avg) <- c("Conc.Air")
}

# Read calculated average sampling rates ----------------------------------
sr <- read.csv("Output/Data/csv/ParticipantSRV02.csv")
# Select only average sampling rate
sr <- sr[, 1:2]

# Select wore WBs ---------------------------------------------------------
{
  wb.Mi.l <- data[4, c(2, 3:175)]
  wb.Mi.r <- data[5, c(2, 3:175)]
  wb.Ea.r <- data[6, c(2, 3:175)]
  wb.Ea.l <- data[7, c(2, 3:175)]
  wb.Ya.r <- data[8, c(2, 3:175)]
  wb.Ya.l <- data[9, c(2, 3:175)]
  wb.An.r <- data[10, c(2, 3:175)]
  wb.An.l <- data[11, c(2, 3:175)]
  wb.Xu.r <- data[12, c(2, 3:175)]
  wb.Xu.l <- data[13, c(2, 3:175)]
}
# Combined wore WBs
wb.wr <- rbind(wb.Mi.l, wb.Mi.r, wb.Ea.r, wb.Ea.l, wb.Ya.r, wb.Ya.l,
               wb.An.r, wb.An.l, wb.Xu.r, wb.Xu.l)

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
# !!! Use a constant sampling rate
# wb_div_sr <- sweep(wb_common, 2, 1.5, FUN = "/")
# Extract time.day from wb.wr
wb_time_day <- wb.wr$time.day
# Divide wb_div_sr further by corresponding value in wb_time_day
conc.wb <- sweep(wb_div_sr, 1, wb_time_day, FUN = "/")
rownames(conc.wb) <- c('wb.Mi.l', 'wb.Mi.r', 'wb.Ea.r', 'wb.Ea.l',
                       'wb.Ya.r', 'wb.Ya.l', 'wb.An.l', 'wb.An.r',
                       'wb.Xu.r', 'wb.Xu.l')

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
  Air_Concentration = rep(tPCB.conc.air, times = c(4, 2, 2, 2, 2,2)),
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
  TRUE ~ NA_character_  # This handles any unmatched cases
)

# Define a color palette with enough distinct colors for the number of volunteers
color_palette <- c(
  "#377eb8", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#e41a1c", "#264c73",
  "#f781bf", "#a65628")

# Define a shape palette with enough distinct shapes for the number of volunteers
shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Create the plot
plotAirWBtPCB <- ggplot(data, aes(x = Air_Concentration, y = Wb_Concentration,
                           fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 3.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(1, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x),
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
ggsave("Output/Plots/AirWBtPCBSR.png", plot = plotAirWBtPCB, width = 6,
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

# Calculate mean percent error
mean_error <- mean(percentage_error)
print(paste("Mean Error:", mean_error))

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
ggsave("Output/Plots/AirWBPCBi.png", plot = plotAirWBPCBi, width = 6,
       height = 5, dpi = 500)


