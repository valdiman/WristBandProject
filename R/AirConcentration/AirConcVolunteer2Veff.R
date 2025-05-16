## Concentration estimation
# Office and home


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

# Read data -----------------------------------------
{
  data.2 <- read.csv("Data/VolunteersV02.csv")
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

# Have the same number of congeners
common_ids <- intersect(ko$congener, logKwb$congener)
ko <- ko[ko$congener %in% common_ids, ]
logKwb <- logKwb[logKwb$congener %in% common_ids, ]
data.2.1 <- data.2[, common_ids]

# Calculate air PCB concentration from static WBs -------------------------
{
  # Calculate effective volume for static WBs
  # 2 offices, 1 and 2 are replicates
  # Use effective volume. Adult WBs
  Vwb.of1.1 <- data.2$vol.WB[1] # [m3]
  Awb.of1.1 <- data.2$area.WB[1] # [m2]
  veff_static.of1.1 <- 10^(logKwb$logKwb) * Vwb.of1.1 * 
    (1 - exp(-ko$ko * Awb.of1.1 / Vwb.of1.1 / 10^(logKwb$logKwb) * data.2$office.time.day[1]))
  Vwb.of1.2 <- data.2$vol.WB[2] # [m3]
  Awb.of1.2 <- data.2$area.WB[2] # [m2]
  veff_static.of1.2 <- 10^(logKwb$logKwb) * Vwb.of1.2 * 
    (1 - exp(-ko$ko * Awb.of1.2 / Vwb.of1.2 / 10^(logKwb$logKwb) * data.2$office.time.day[2]))
  Vwb.of2 <- data.2$vol.WB[3] # [m3]
  Awb.of2 <- data.2$area.WB[3] # [m2]
  veff_static.of2 <- 10^(logKwb$logKwb) * Vwb.of2 * 
    (1 - exp(-ko$ko * Awb.of2 / Vwb.of2 / 10^(logKwb$logKwb) * data.2$office.time.day[3]))

  # Calculate air concentration in ng/m3 from static WBs
  # For of1.1
  conc.of1.1 <- as.data.frame(t(data.2.1[1, ] / veff_static.of1.1))
  # For of1.2
  conc.of1.2 <- as.data.frame(t(data.2.1[2, ] / veff_static.of1.2))
  # Bind the data frames together row-wise
  conc.of1 <- cbind(conc.of1.1, conc.of1.2)
  # Average 1.1 & 1.2
  conc.of1 <- as.data.frame(colMeans(t(conc.of1)))
  # Bind the data frames together row-wise
  # For of2
  conc.of2 <- as.data.frame(t(data.2.1[3, ] / veff_static.of2))
  # Combine concentrations
  conc.air <- cbind(conc.of1, conc.of2)
  # Change column names of the last three columns
  colnames(conc.air) <- c("Conc.Air.1", "Conc.Air.2")
}

# Calculate air concentration from worn WBs -------------------------------
# Read ko from PersonalSamplingRatesV02
ko.p <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSRV02.csv")
ko.p <- ko.p[c(1,7)]

{
  Vwb.V1 <- data.2$vol.WB[4] # [m3]
  Awb.V1 <- data.2$area.WB[4] # [m2]
  veff.V1.o <- 10^(logKwb$logKwb) * Vwb.V1 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V1 / Vwb.V1 / 10^(logKwb$logKwb) * data.2$office.time.day[4]))
  veff.V1.h <- 10^(logKwb$logKwb) * Vwb.V1 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V1 / Vwb.V1 / 10^(logKwb$logKwb) * data.2$total.time.day[5]))
  Vwb.V2 <- data.2$vol.WB[6] # [m3]
  Awb.V2 <- data.2$area.WB[6] # [m2]
  veff.V2.o <- 10^(logKwb$logKwb) * Vwb.V2 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V2 / Vwb.V2 / 10^(logKwb$logKwb) * data.2$office.time.day[7]))
  veff.V2.h <- 10^(logKwb$logKwb) * Vwb.V2 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V2 / Vwb.V2 / 10^(logKwb$logKwb) * data.2$total.time.day[6]))
  Vwb.V3 <- data.2$vol.WB[8] # [m3]
  Awb.V3 <- data.2$area.WB[8] # [m2]
  veff.V3.o <- 10^(logKwb$logKwb) * Vwb.V3 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V3 / Vwb.V3 / 10^(logKwb$logKwb) * data.2$office.time.day[8]))
  veff.V3.h <- 10^(logKwb$logKwb) * Vwb.V3 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V3 / Vwb.V3 / 10^(logKwb$logKwb) * data.2$total.time.day[9]))
  Vwb.V4 <- data.2$vol.WB[10] # [m3]
  Awb.V4 <- data.2$area.WB[10] # [m2]
  veff.V4.o <- 10^(logKwb$logKwb) * Vwb.V4 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4 / Vwb.V4 / 10^(logKwb$logKwb) * data.2$office.time.day[11]))
  veff.V4.h <- 10^(logKwb$logKwb) * Vwb.V4 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V4 / Vwb.V4 / 10^(logKwb$logKwb) * data.2$total.time.day[10]))
  Vwb.V5 <- data.2$vol.WB[12] # [m3]
  Awb.V5 <- data.2$area.WB[12] # [m2]
  veff.V5.o <- 10^(logKwb$logKwb) * Vwb.V5 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5 / Vwb.V5 / 10^(logKwb$logKwb) * data.2$office.time.day[13]))
  veff.V5.h <- 10^(logKwb$logKwb) * Vwb.V5 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V5 / Vwb.V5 / 10^(logKwb$logKwb) * data.2$total.time.day[12]))
}

# Estimate air concentration in ng/m3 from WBs
# Using Veff
{
  conc.V1.o <- as.data.frame(t(data.2.1[4, ] / veff.V1.o))
  conc.V1.h <- as.data.frame(t(data.2.1[5, ] / veff.V1.h))
  conc.V2.o <- as.data.frame(t(data.2.1[7, ] / veff.V2.o))
  conc.V2.h <- as.data.frame(t(data.2.1[6, ] / veff.V2.h))
  conc.V3.o <- as.data.frame(t(data.2.1[8, ] / veff.V3.o))
  conc.V3.h <- as.data.frame(t(data.2.1[9, ] / veff.V3.h))
  conc.V4.o <- as.data.frame(t(data.2.1[11, ] / veff.V4.o))
  conc.V4.h <- as.data.frame(t(data.2.1[10, ] / veff.V4.h))
  conc.V5.o <- as.data.frame(t(data.2.1[13, ] / veff.V5.o))
  conc.V5.h <- as.data.frame(t(data.2.1[12, ] / veff.V5.h))
}

# Combined wore WBs
conc.wb <- cbind(conc.V1.o, conc.V1.h, conc.V2.o, conc.V2.h, conc.V3.o,
                 conc.V3.h, conc.V4.o, conc.V4.h, conc.V5.o, conc.V5.h)

colnames(conc.wb) <- c('conc.V1.o', 'conc.V1.h', 'conc.V2.o', 'conc.V2.h',
                       'conc.V3.o', 'conc.V3.h', 'conc.V4.o', 'conc.V4.h',
                       'conc.V5.o', 'conc.V5.h')

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- colSums(conc.wb, na.rm = TRUE)
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Export results
write.csv(conc.air,
          file = "Output/Data/csv/Volunteer/VolunteerConcStaticWB2.csv")
write.csv(conc.wb,
          file = "Output/Data/csv/Volunteer/VolunteerConcWB2.csv")


# Plots -------------------------------------------------------------------
# tPCB
# Create a data frame with the combined data
data.plot <- data.frame(
  Wb_Concentration = tPCB.conc.wb,
  Air_Concentration = rep(tPCB.conc.air, times = c(6, 4)),
  Volunteer = names(tPCB.conc.wb)
)

# Change names for legend
data.plot$Volunteer2 <- case_when(
  data.plot$Volunteer == "conc.V1.o" ~ "Vol. 1 Of",
  data.plot$Volunteer == "conc.V1.h" ~ "Vol. 1 Ho",
  data.plot$Volunteer == "conc.V3.o" ~ "Vol. 2 Of",
  data.plot$Volunteer == "conc.V3.h" ~ "Vol. 2 Ho",
  data.plot$Volunteer == "conc.V2.o" ~ "Vol. 3 Of",
  data.plot$Volunteer == "conc.V2.h" ~ "Vol. 3 Ho",
  data.plot$Volunteer == "conc.V5.o" ~ "Vol. 8 Of",
  data.plot$Volunteer == "conc.V5.h" ~ "Vol. 8 Ho",
  data.plot$Volunteer == "conc.V4.o" ~ "Vol. 9 Of",
  data.plot$Volunteer == "conc.V4.h" ~ "Vol. 9 Ho",
  TRUE ~ NA_character_
)

# Create a new column that will group the data into Vol. 1, Vol. 2, etc.
data.plot$Volunteer_Group <- case_when(
  grepl("conc.V1", data.plot$Volunteer) ~ "Vol. 1",
  grepl("conc.V2", data.plot$Volunteer) ~ "Vol. 3",
  grepl("conc.V3", data.plot$Volunteer) ~ "Vol. 2",
  grepl("conc.V4", data.plot$Volunteer) ~ "Vol. 9",
  grepl("conc.V5", data.plot$Volunteer) ~ "Vol. 8",
  TRUE ~ NA_character_
)

# Reshape the data into long format so that each volunteer will have three rows: Of, Ho, and Air
data_long <- data.plot %>%
  # Create a new "Category" column to represent the Concentration type
  mutate(Concentration_Type = case_when(
    grepl("Of", data.plot$Volunteer2) ~ "Of",
    grepl("Ho", data.plot$Volunteer2) ~ "Ho",
    TRUE ~ "Air"
  )) %>%
  # Reshape to long format
  pivot_longer(cols = c(Wb_Concentration, Air_Concentration), 
               names_to = "Concentration_Name", 
               values_to = "Concentration") %>%
  # Now, add a column to distinguish between the three types of concentration
  mutate(Concentration_Category = case_when(
    Concentration_Name == "Wb_Concentration" & Concentration_Type == "Of" ~ "Of",
    Concentration_Name == "Wb_Concentration" & Concentration_Type == "Ho" ~ "Ho",
    Concentration_Name == "Air_Concentration" ~ "Air",
    TRUE ~ NA_character_
  ))

# Set the factor levels for Concentration_Category to control the order
data_long$Concentration_Category <- factor(
  data_long$Concentration_Category,
  levels = c("Air", "Of", "Ho")  # Desired order
)

# Plot
plotAirWBtPCB <- ggplot(data_long, aes(x = Volunteer_Group, y = Concentration,
                      fill = Concentration_Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Dodge for side-by-side bars
  scale_fill_manual(
    values = c("Ho" = "#E69F00", "Of" = "blue", "Air" = "#009E73"),  # Custom colors
    labels = c("Ho" = "WB full-day", "Of" = "WB office-only",
               "Air" = "Air PCB office")  # Custom labels
  ) +
  labs(x = '', fill = "Location") +
  ylab(expression(bold("Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14,
                                   angle = 45, hjust =1),
        axis.title.x = element_text(face = "bold", size = 14))

# Print the plot
print(plotAirWBtPCB)

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/AirWBtPCBOfficeHomeVeff.png",
       plot = plotAirWBtPCB, width = 6, height = 5, dpi = 500)

# Plot 1:1 Air PCB office vs. WB full-day and WB only office

# Define a color palette with enough distinct colors for the number of volunteers
color_palette <- c("#377eb8", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#f781bf", "#a65628", "#999999", "#e41a1c", "#4daf4a")

# Define a shape palette with enough distinct shapes for the number of volunteers
shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25)

# Create the plot
plot <- ggplot(data.plot, aes(x = Air_Concentration, y = Wb_Concentration,
                              fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 3.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(1, 10^2),
                breaks = 10^(0:2),  # Integer powers of 10 only
                labels = scales::trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^2),
                breaks = 10^(0:2),  # Integer powers of 10 only
                labels = scales::trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Office Air Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  guides(fill = guide_legend(override.aes = list(color = "black", shape = shape_palette))) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(fill = "Volunteers", shape = "Volunteers") +
  theme(legend.position = "right")

# See plot
plot

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VolunteerAirOffHomeVeff.png",
       plot = plot, width = 6, height = 5, dpi = 500)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
data.plot$factor2 <- data.plot$Wb_Concentration/data.plot$Air_Concentration

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(data.plot[data.plot$factor2 > 0.5 & data.plot$factor2 < 2, ])/nrow(data.plot)*100

# Estimate percentage error ---------------------------------------------------
# Calculate percentage error
percentage_error <-function(observed, predicted) {
  return(abs(observed - predicted)/abs(observed) * 100)
}

# Calculate percentage errors
percentage_errors <- percentage_error(data.plot$Air_Concentration,
                                     data.plot$Wb_Concentration)

# Calculate mean percent error
mean_error <- mean(percentage_errors)
print(paste("Mean Error:", mean_error))

# Ratios between home & office ------------------------------------------
# Select data
wb_office <- tPCB.conc.wb[grepl("\\.o$", names(tPCB.conc.wb))]
wb_home <- tPCB.conc.wb[grepl("\\.h$", names(tPCB.conc.wb))] 
# Remove .o and .h from names
names(wb_office) <- gsub("\\.o$", "", names(wb_office))
names(wb_home) <- gsub("\\.h$", "", names(wb_home))
# Match the names to ensure alignment
common_volunteers <- intersect(names(wb_office), names(wb_home))
# Subset the vectors to include only common volunteers
wb_office <- wb_office[common_volunteers]
wb_home <- wb_home[common_volunteers]
# Calculate the ratio of Office to Home concentrations
ratios <- wb_office / wb_home
# Create a data frame for comparison including the ratios
comparison_df <- data.frame(
  Volunteer = common_volunteers,
  WB_Office = wb_office,
  WB_Home = wb_home,
  Ratio = ratios,
  stringsAsFactors = FALSE
)
# Print the result
print(comparison_df)

# Plot individual PCB congeners -------------------------------------------
# Air Add PCB as a column
conc_air_common <- conc.air %>% rownames_to_column(var = "PCB")

# Transform conc_air_common into long format
conc_air_long <- conc_air_common %>%
  pivot_longer(
    cols = starts_with("Conc.Air"),  # This will select Conc.Air.1 and Conc.Air.2
    names_to = "Conc.Type",
    values_to = "Conc.Air"
  )

# Replicate the rows: 3 times for Conc.Air.1 and 2 times for Conc.Air.2
conc_air_long <- conc_air_long %>%
  mutate(Replicate_Count = case_when(
    Conc.Type == "Conc.Air.1" ~ 3,  # Replicate Conc.Air.1 3 times
    Conc.Type == "Conc.Air.2" ~ 2,  # Replicate Conc.Air.2 2 times
    TRUE ~ 1
  )) %>%
  uncount(weights = Replicate_Count) %>%  # Replicate based on the count
  mutate(Volunteer = case_when(
    Conc.Type == "Conc.Air.1" ~ rep(c("Vol1", "Vol2", "Vol3"), length.out = n()),  # Assign Vol1, Vol2, Vol3
    Conc.Type == "Conc.Air.2" ~ rep(c("Vol4", "Vol5"), length.out = n()),  # Assign Vol4, Vol5
    TRUE ~ NA_character_
  ))

# Transpose conc.wb
t.conc.wb <- t(conc.wb)

# WB Move row names to the first column
conc_wb_clean <- t.conc.wb %>%
  as.data.frame() %>%
  rownames_to_column(var = "Row.Name")  # Row names are moved to 'Row.Name' column

# Remove rows that end with 'o'
conc_wb_clean <- conc_wb_clean %>%
  filter(!grepl("o$", Row.Name))  # Filter out rows where the 'Row.Name' ends with 'o'

# Rename the first column to 'Volunteer'
conc_wb_clean <- conc_wb_clean %>%
  rename(Volunteer = Row.Name)

# Assign 'Vol1', 'Vol2', ..., 'Vol5' based on row numbers
conc_wb_clean <- conc_wb_clean %>%
  mutate(Volunteer = paste0("Vol", rep(1:5, length.out = n())))  # Assigns Vol1 to Vol5 cyclically

conc_wb_clean_long <- conc_wb_clean %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "PCB", 
               values_to = "Conc.WB")

# Merge the two data frames
merged_data <- conc_air_long %>%
  inner_join(conc_wb_clean_long, by = c("PCB", "Volunteer" = "Volunteer"))

# Change column name from PCB to congener
merged_data <- merged_data %>%
  rename(congener = PCB)

# Filter out rows with 0 in either Conc.Air or Conc.WB
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

# Add Ho after volunteer
filtered_data$Volunteer <- paste0(filtered_data$Volunteer, " Ho")

# Color and shapes from tPCB plot
# Define a color palette with enough distinct colors for the number of volunteers
color_palette <- c("#377eb8", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

# Define a shape palette with enough distinct shapes for the number of volunteers
shape_palette <- c(21, 22, 23, 24, 25)

# Plot with different colors for each volunteer  
p.AirWBPCBi.volun2 <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
                                           fill = Volunteer, shape = Volunteer)) +
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
  xlab(expression(bold("Office Air Concentration PCBi (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted full-day Concentration PCBi (ng/m"^3*")"))) +
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
  theme(legend.position = "right")

# See plot
p.AirWBPCBi.volun2

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VoluntHomeOfficeAirWBPCBVeff.png",
       plot = p.AirWBPCBi.volun2, width = 6, height = 5, dpi = 500)

# Prepare data to combine with Frederiksen
Volunteer2_PCBi_fred <- filtered_data %>%
  filter(congener %in% c("PCB8", "PCB18.30", "PCB20.28", "PCB31", "PCB44.47.65",
                    "PCB52", "PCB66", "PCB99", "PCB90.101.113", "PCB105",
                    "PCB118", "PCB129.138.163", "PCB153.168", "PCB180.193"))

Volunteer2_PCBi_fred <- Volunteer2_PCBi_fred %>%
  mutate(congener = factor(congener, levels = unique(congener)))

# Plot
color_palette2 <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette2 <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
p.AirWBPCBi.volun2.Fred <- ggplot(Volunteer2_PCBi_fred, aes(x = Conc.Air, y = Conc.WB, 
                                           fill = congener, shape = congener)) +
  geom_point(size = 2.5, color = "black", stroke = 0.5) + # Points will have black border
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(0.001, 10^4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Office Air Concentration PCBi (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted Concentration PCBi (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  scale_fill_manual(values = color_palette2) +
  scale_shape_manual(values = shape_palette2) +
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

# See plot
p.AirWBPCBi.volun2.Fred

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/VoluntHomeOfficeAirWBPCBFredVeff.png",
       plot = p.AirWBPCBi.volun2.Fred, width = 6, height = 5, dpi = 500)

# Export data
write.csv(Volunteer2_PCBi_fred,
          file = "Output/Data/csv/FrederiksenPCB/Volunteer2_PCBiVeff.csv")


