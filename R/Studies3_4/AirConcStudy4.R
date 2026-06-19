## Concentration estimation for Study 4
# Office and home


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

# Read data -----------------------------------------
{
  data.s4 <- read.csv("Data/IRO/11_SampleWBMassStudy4.csv", check.names = FALSE)
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

# Select only congeners
common_ids <- intersect(ko$congener, logKwb$congener)
ko <- ko[ko$congener %in% common_ids, ]
logKwb <- logKwb[logKwb$congener %in% common_ids, ]
data.s4.1 <- data.s4[, common_ids]

# Calculate air PCB concentration from static WBs -------------------------
{
  # Calculate effective volume for static WBs
  # V15, 16 and 17 in the same office (2 replicates, but different WBs type)
  Vwb.1.V15.16.17 <- data.s4$vol.WB[1] # [m3]
  Awb.1.V15.16.17 <- data.s4$area.WB[1] # [m2]
  veff_static.1.V15.16.17 <- 10^(logKwb$logKwb) * Vwb.1.V15.16.17 * 
    (1 - exp(-ko$ko * Awb.1.V15.16.17 / Vwb.1.V15.16.17 / 10^(logKwb$logKwb) * data.s4$time[1] / 24))
  Vwb.2.V15.16.17 <- data.s4$vol.WB[2] # [m3]
  Awb.2.V15.16.17 <- data.s4$area.WB[2] # [m2]
  veff_static.2.V15.16.17 <- 10^(logKwb$logKwb) * Vwb.2.V15.16.17 * 
    (1 - exp(-ko$ko * Awb.2.V15.16.17 / Vwb.2.V15.16.17 / 10^(logKwb$logKwb) * data.s4$time[2] / 24))
  # Volunteers 18 and 19 same office
  Vwb.V18.19 <- data.s4$vol.WB[3] # [m3]
  Awb.V18.19 <- data.s4$area.WB[3] # [m2]
  veff_static.V18.19 <- 10^(logKwb$logKwb) * Vwb.V18.19 * 
    (1 - exp(-ko$ko * Awb.V18.19 / Vwb.V18.19 / 10^(logKwb$logKwb) * data.s4$time[3] / 24))

  # Calculate air concentration in ng/m3 from static WBs
  # For V15.16.17 (1)
  conc.1.V15.16.17 <- as.data.frame(t(data.s4.1[1, ] / veff_static.1.V15.16.17))
  # For V15.16.17 (2)
  conc.2.V15.16.17 <- as.data.frame(t(data.s4.1[2, ] / veff_static.2.V15.16.17))
  # For V15.16.17
  conc.V15.16.17 <- (conc.1.V15.16.17 + conc.2.V15.16.17) / 2
  # For V15
  conc.V15 <- conc.V15.16.17
  # For V16
  conc.V16 <- conc.V15.16.17
  # For V17
  conc.V17 <- conc.V15.16.17
  # For V18
  conc.V18 <- as.data.frame(t(data.s4.1[3, ] / veff_static.V18.19))
  # For V19
  conc.V19 <- as.data.frame(t(data.s4.1[3, ] / veff_static.V18.19))
  # Combine concentrations
  conc.air <- cbind(conc.V15, conc.V16, conc.V17, conc.V18, conc.V19)
  # Change column names of the last three columns
  colnames(conc.air) <- c("Conc.Air.V15", "Conc.Air.V16", "Conc.Air.V17",
                          "Conc.Air.V18", "Conc.Air.V19")
}

# Check total PCB
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)
# See
tPCB.conc.air

# Calculate air concentration from worn WBs -------------------------------
# Read ko from PersonalSamplingRates.R
ko.p <- read.csv("Output/Data/Studies1_2/Study2/PersonalAveSR.csv")
ko.p <- ko.p[c(1,7)]
{
  Vwb.V15 <- data.s4$vol.WB[4] # [m3]
  Awb.V15 <- data.s4$area.WB[4] # [m2]
  veff.V15.o <- 10^(logKwb$logKwb) * Vwb.V15 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V15 / Vwb.V15 / 10^(logKwb$logKwb) * data.s4$time[4] / 24))
  veff.V15.d <- 10^(logKwb$logKwb) * Vwb.V15 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V15 / Vwb.V15 / 10^(logKwb$logKwb) * data.s4$time[5] / 24))
  Vwb.V16 <- data.s4$vol.WB[6] # [m3]
  Awb.V16 <- data.s4$area.WB[6] # [m2]
  veff.V16.o <- 10^(logKwb$logKwb) * Vwb.V16 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V16 / Vwb.V16 / 10^(logKwb$logKwb) * data.s4$time[6] / 24))
  veff.V16.d <- 10^(logKwb$logKwb) * Vwb.V16 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V16 / Vwb.V16 / 10^(logKwb$logKwb) * data.s4$time[7] / 24))
  Vwb.V17 <- data.s4$vol.WB[8] # [m3]
  Awb.V17 <- data.s4$area.WB[8] # [m2]
  veff.V17.o <- 10^(logKwb$logKwb) * Vwb.V17 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V17 / Vwb.V17 / 10^(logKwb$logKwb) * data.s4$time[8] / 24))
  veff.V17.d <- 10^(logKwb$logKwb) * Vwb.V17 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V17 / Vwb.V17 / 10^(logKwb$logKwb) * data.s4$time[9] / 24))
  Vwb.V18 <- data.s4$vol.WB[10] # [m3]
  Awb.V18 <- data.s4$area.WB[10] # [m2]
  veff.V18.o <- 10^(logKwb$logKwb) * Vwb.V18 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V18 / Vwb.V18 / 10^(logKwb$logKwb) * data.s4$time[10] / 24))
  veff.V18.d <- 10^(logKwb$logKwb) * Vwb.V18 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V18 / Vwb.V18 / 10^(logKwb$logKwb) * data.s4$time[11] / 24))
  Vwb.V19 <- data.s4$vol.WB[12] # [m3]
  Awb.V19 <- data.s4$area.WB[12] # [m2]
  veff.V19.o <- 10^(logKwb$logKwb) * Vwb.V19 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V19 / Vwb.V19 / 10^(logKwb$logKwb) * data.s4$time[12] / 24))
  veff.V19.d <- 10^(logKwb$logKwb) * Vwb.V19 * 
    (1 - exp(-ko.p$Average_ko2 * Awb.V19 / Vwb.V19 / 10^(logKwb$logKwb) * data.s4$time[13] / 24))
}

# Estimate air concentration in ng/m3 from WBs
# Using Veff
{
  conc.V15.o <- as.data.frame(t(data.s4.1[4, ] / veff.V15.o))
  conc.V15.d <- as.data.frame(t(data.s4.1[5, ] / veff.V15.d))
  conc.V16.o <- as.data.frame(t(data.s4.1[6, ] / veff.V16.o))
  conc.V16.d <- as.data.frame(t(data.s4.1[7, ] / veff.V16.d))
  conc.V17.o <- as.data.frame(t(data.s4.1[8, ] / veff.V17.o))
  conc.V17.d <- as.data.frame(t(data.s4.1[9, ] / veff.V17.d))
  conc.V18.o <- as.data.frame(t(data.s4.1[10, ] / veff.V18.o))
  conc.V18.d <- as.data.frame(t(data.s4.1[11, ] / veff.V18.d))
  conc.V19.o <- as.data.frame(t(data.s4.1[12, ] / veff.V19.o))
  conc.V19.d <- as.data.frame(t(data.s4.1[13, ] / veff.V19.d))
}

# Combined wore WBs
conc.wb <- cbind(conc.V15.o, conc.V15.d, conc.V16.o, conc.V16.d, conc.V17.o,
                 conc.V17.d, conc.V18.o, conc.V18.d, conc.V19.o, conc.V19.d)
colnames(conc.wb) <- c('Conc.WB.V15.o', 'Conc.WB.V15.d', 'Conc.WB.V16.o',
                       'Conc.WB.V16.d', 'Conc.WB.V17.o', 'Conc.WB.V17.d',
                       'Conc.WB.V18.o', 'Conc.WB.V18.d', 'Conc.WB.V19.o',
                       'Conc.WB.V19.d')

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- colSums(conc.wb, na.rm = TRUE)
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Export results
write.csv(conc.air,
          file = "Output/Data/Study4/VolunteerConcStaticWBStudy4.csv")
write.csv(conc.wb,
          file = "Output/Data/Study4/VolunteerConcWBStudy4.csv")

# Plots -------------------------------------------------------------------
# tPCB
# Create a data frame with the combined data
data.conc <- data.frame(
  Wb_Concentration = tPCB.conc.wb,
  Air_Concentration = rep(tPCB.conc.air, each = 2),
  Volunteer = names(tPCB.conc.wb)
)

# Change names for legend
data.conc$Volunteer2 <- case_when(
  data.conc$Volunteer == "Conc.WB.V15.o" ~ "V15 off",
  data.conc$Volunteer == "Conc.WB.V15.d" ~ "V15 day",
  data.conc$Volunteer == "Conc.WB.V16.o" ~ "V16 off",
  data.conc$Volunteer == "Conc.WB.V16.d" ~ "V16 day",
  data.conc$Volunteer == "Conc.WB.V17.o" ~ "V17 off",
  data.conc$Volunteer == "Conc.WB.V17.d" ~ "V17 day",
  data.conc$Volunteer == "Conc.WB.V18.o" ~ "V18 off",
  data.conc$Volunteer == "Conc.WB.V18.d" ~ "V18 day",
  data.conc$Volunteer == "Conc.WB.V19.o" ~ "V19 off",
  data.conc$Volunteer == "Conc.WB.V19.d" ~ "V19 day",
  TRUE ~ NA_character_
)

# To organize them
data.conc$Volunteer2 <- factor(
  data.conc$Volunteer2,
  levels = str_sort(unique(data.conc$Volunteer2), numeric = TRUE))

color_palette <- c("#377eb8", "#377eb8", "#ff7f0e", "#ff7f0e", "#2ca02c",
                   "#2ca02c", "#d62728", "#d62728", "#9467bd", "#9467bd")

shape_palette <- c(21, 24, 21, 24, 21, 24, 21, 24, 21, 24)

# Create the plot
plot.s4 <- ggplot(data.conc, aes(x = Air_Concentration, y = Wb_Concentration,
                                       fill = Volunteer2, shape = Volunteer2)) +
  geom_point(size = 3.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(1, 10^2), breaks = 10^(0:3),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^2), breaks = 10^(0:3),
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
print(plot.s4)

# Save plot in folder
ggsave("Output/Plots/Study4/VolunteerAirOffHomeVeffStudy4.png",
       plot = plot.s4, width = 6, height = 5, dpi = 500)

# Bar plot
# Format data
data_bar <- data.conc %>%
  pivot_longer(
    cols = c(Wb_Concentration, Air_Concentration),
    names_to = "Concentration_Name", values_to = "Concentration"
  ) %>%
  mutate(
    Concentration_Category = case_when(
      Concentration_Name == "Air_Concentration" ~ "Air",
      grepl("off", Volunteer2) ~ "off",
      grepl("day", Volunteer2) ~ "day"),
    Volunteer_Group = sub(" (off|day)$", "", Volunteer2))

# Plot
plot.s4.2 <- ggplot(data_bar, aes(x = Volunteer_Group, y = Concentration,
                                  fill = Concentration_Category)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black") +
  scale_fill_manual(
    values = c("day" = "#E69F00",
               "off" = "#0072B2",
               "Air" = "#009E73"),
    labels = c("day" = "WB full-day",
               "off" = "WB office-only",
               "Air" = "Air PCB office")) +
  labs(x = "", fill = "Measurement",
       y = expression(bold("Concentration " * Sigma * "PCB (ng/m"^3*")"))) +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14,
                               angle = 45, hjust = 1),
    axis.title.x = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12))

# Print the plot
print(plot.s4.2)

# Save plot in folder
ggsave("Output/Plots/Study4/BarVolunteerAirOffHomeVeffStudy4.png",
       plot = plot.s4.2, width = 6, height = 5, dpi = 500)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
data.conc$factor2 <- signif(data.conc$Wb_Concentration/data.conc$Air_Concentration, 2)

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(data.conc[data.conc$factor2 >= 0.5 & data.conc$factor2 <= 2, ])/nrow(data.conc)*100

# Estimate percentage error ---------------------------------------------------
# Calculate percentage error
percentage_error <-function(observed, predicted) {
  return(abs(observed - predicted)/abs(observed) * 100)
}

# Calculate percentage errors
percentage_errors <- percentage_error(data.conc$Air_Concentration,
                                     data.conc$Wb_Concentration)

# Calculate mean percent error
mean_error <- mean(percentage_errors)
print(paste("Mean Error:", mean_error))

# Ratios between home & office ------------------------------------------
# Select data
wb_office <- tPCB.conc.wb[grepl("\\.o$", names(tPCB.conc.wb))]
wb_day <- tPCB.conc.wb[grepl("\\.d$", names(tPCB.conc.wb))] 
# Remove .o and .h from names
names(wb_office) <- gsub("\\.o$", "", names(wb_office))
names(wb_day) <- gsub("\\.d$", "", names(wb_day))
# Match the names to ensure alignment
common_volunteers <- intersect(names(wb_office), names(wb_day))
# Subset the vectors to include only common volunteers
wb_office <- wb_office[common_volunteers]
wb_day <- wb_day[common_volunteers]
# Calculate the ratio of Office to Home concentrations
ratios <- wb_office / wb_day
# Create a data frame for comparison including the ratios
comparison_df <- data.frame(
  Volunteer = common_volunteers,
  WB_Office = wb_office,
  WB_Day = wb_day,
  Ratio = ratios,
  stringsAsFactors = FALSE
)
# Print the result
print(comparison_df)

# Plot individual PCB congeners -------------------------------------------
# This is just for visualization purposes
# Create a data frame with the combined data
conc_air_long <- conc.air %>%
  rownames_to_column("PCB") %>%
  pivot_longer(-PCB, names_to = "Volunteer", values_to = "Conc.Air") %>%
  mutate(Volunteer = str_extract(Volunteer, "V[0-9]+"))

conc_wb_long <- conc.wb %>%
  rownames_to_column("PCB") %>%
  pivot_longer(-PCB, names_to = "Volunteer_full", values_to = "Conc.WB") %>%
  mutate(
    Volunteer = str_extract(Volunteer_full, "V[0-9]+"),
    Condition = ifelse(grepl("\\.o$", Volunteer_full), "Office WB", "Full-day WB"))

# Merge
filtered_data <- conc_air_long %>%
  inner_join(conc_wb_long, by = c("PCB", "Volunteer")) %>%
  filter(Conc.Air > 0, Conc.WB > 0)

# Factor
filtered_data$Volunteer <- factor(filtered_data$Volunteer,
                                  levels = paste0("V", 15:19),
                                  labels = paste0("V", 15:19))

# Colors
color_palette <- c("#377eb8", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
shape_palette <- c(21, 22, 23, 24, 25)

# Plot
plot.s4.3 <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB,
                          fill = Volunteer, shape = Volunteer)) +
  geom_point(size = 2.5, color = "black", stroke = 0.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  scale_y_log10(
    limits = c(1e-5, 1e2),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(
    limits = c(1e-5, 1e2),
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

# See plot
print(plot.s4.3)

# Save plot in folder
ggsave("Output/Plots/Study4/VoluntHomeOfficeAirWBPCBVeffStudy4.png",
       plot = plot.s4.3, width = 8, height = 5, dpi = 500)

# Estimate error (factor) --------------------------------------------
# below 1:1 → model underestimates WB
# above 1:1 → model overestimates WB
# within factor 2 → acceptable agreement

factor_summary <- filtered_data %>%
  mutate(factor2 = Conc.WB / Conc.Air) %>%
  group_by(Condition) %>%
  summarise(
    n = n(),
    pct_within_2 = mean(factor2 >= 0.5 & factor2 <= 2, na.rm = TRUE) * 100,
    pct_below_1to1 = mean(factor2 < 1, na.rm = TRUE) * 100,
    pct_above_1to1 = mean(factor2 > 1, na.rm = TRUE) * 100)

factor_summary

# Prepare data to combine with Frederiksen
Volunteer2_PCBi_fred <- filtered_data %>%
  filter(PCB %in% c("PCB8", "PCB18+30", "PCB20+28", "PCB31", "PCB44+47+65",
                    "PCB52", "PCB66", "PCB99", "PCB90+101+113", "PCB105",
                    "PCB118", "PCB129+138+163", "PCB153+168", "PCB180+193"))

Volunteer2_PCBi_fred <- Volunteer2_PCBi_fred %>%
  mutate(PCB = factor(PCB, levels = unique(PCB)))

# Plot
color_palette2 <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette2 <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
p.AirWBPCBi.volun2.Fred <- ggplot(Volunteer2_PCBi_fred,
       aes(x = Conc.Air, y = Conc.WB,
           fill = PCB, shape = PCB)) +
  geom_point(size = 2.5, color = "black", stroke = 0.4, alpha = 0.7) +
  scale_fill_manual(values = color_palette2) +
  scale_shape_manual(values = shape_palette2) +
  scale_y_log10(
    limits = c(1e-4, 1e2),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(
    limits = c(1e-4, 1e2),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
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
  labs(fill = "PCB", shape = "PCB")

# See plot
p.AirWBPCBi.volun2.Fred

# Save plot in folder
ggsave("Output/Plots/Study4/VoluntHomeOfficeAirWBPCBFredCongenersVeffStudy4.png",
       plot = p.AirWBPCBi.volun2.Fred, width = 8, height = 5, dpi = 500)

# Export data
write.csv(Volunteer2_PCBi_fred,
          file = "Output/Data/Frederiksen/Study4_PCBiVeff.csv")

