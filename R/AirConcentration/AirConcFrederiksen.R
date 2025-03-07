## Script to compare mass accumulated in WBs and air concentration
# in apartments from Frederiksen et al 2022

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("RColorBrewer")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
  library(scales)
}

# Read data from excel ----------------------------------------------------
sr <- data.frame(read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSRV01.csv"))
data.Fred <- data.frame(read_excel("Data/Frederiksen.xlsx", sheet = "Sheet2",
                             col_names = TRUE, col_types = NULL))

# Calculate predicted concentration from WB mass --------------------------
congener_names <- colnames(data.Fred)[6:ncol(data.Fred)]
# Filter sr to keep only matching congeners
sr_filtered <- sr[sr$congener %in% congener_names, ]
# Ensure the order of sr_filtered matches the order in data.Fred
sr_filtered <- sr_filtered[match(congener_names, sr_filtered$congener), ]
# Extract only the Average_Sampling_Rate values
sr.fred <- sr_filtered$Average_Sampling_Rate
# Name the vector with congener names for reference
names(sr.fred) <- sr_filtered$congener
# Select WB masses
wb.Fred <- data.Fred %>% filter(measurement == "mass")
# Find matching columns
matching_cols <- intersect(names(wb.Fred)[6:ncol(wb.Fred)], names(sr.fred))
# Select only relevant columns from wb.Fred for wb.Fred.7d
wb.Fred.7d <- wb.Fred[, c("ID", "group", matching_cols)]
for (col in matching_cols) {
  wb.Fred.7d[[col]] <- wb.Fred[[col]] / (sr.fred[col] * 7)
}
wb.Fred.home <- wb.Fred[, c("ID", "group", matching_cols)]
for (col in matching_cols) {
  wb.Fred.home[[col]] <- wb.Fred[[col]] / (sr.fred[col] * wb.Fred$time.home / 24)
}

# Plot
# Format data
conc.Fred <- data.Fred %>% filter(measurement == "concentration")
# Pivot conc.Fred to long format
conc.Fred_long <- conc.Fred %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.Air")

# 7 days. Pivot wb.Fred.7d to long format
wb.Fred.7d_long <- wb.Fred.7d %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.WB")

# Merge both data frames by the congener column
merged_data <- merge(conc.Fred_long, wb.Fred.7d_long, by = c("ID", "congener"))

# Filter out rows where Conc.Air or Conc.WB are zero
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

filtered_data <- filtered_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Add sampling rates
sr_labels <- paste0(names(sr.fred), " (", round(sr.fred, 2), ")")

# Plot the data
plotfred7d <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
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
  scale_fill_manual(values = color_palette, labels = sr_labels) +
  scale_shape_manual(values = shape_palette, labels = sr_labels) +
  labs(fill = "Congener (Sampling Rate)", shape = "Congener (Sampling Rate)") +
  theme(legend.position = "right")

plotfred7d

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen7d.png", plot = plotfred7d,
       width = 7, height = 7, dpi = 500)

# Home time. Pivot wb.Fred.7d to long format
wb.Fred.home_long <- wb.Fred.home %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.WB")

# Merge both data frames by the congener column
merged_data <- merge(conc.Fred_long, wb.Fred.home_long, by = c("ID", "congener"))

# Filter out rows where Conc.Air or Conc.WB are zero
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

filtered_data <- filtered_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
plotfredhome <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
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
  scale_fill_manual(values = color_palette, labels = sr_labels) +
  scale_shape_manual(values = shape_palette, labels = sr_labels) +
  labs(fill = "Congener (Sampling Rate)", shape = "Congener (Sampling Rate)") +
  theme(legend.position = "right")

plotfredhome

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksenhome.png", plot = plotfredhome,
       width = 7, height = 7, dpi = 500)

# Including sleeping time -------------------------------------------------
# Define the number of sleeping days
sleep_days <- 2  # 8 hrs per day, 6 days
wb.Fred.7d.sl <- wb.Fred[, c("ID", "group", matching_cols)]
# Loop through matching columns
for (col in matching_cols) {
  wb.Fred.7d[[col]] <- wb.Fred[[col]] / (sleep_days * 0.5 + (7 - sleep_days) * sr.fred[col])
}
wb.Fred.home.sl <- wb.Fred[, c("ID", "group", matching_cols)]
for (col in matching_cols) {
  wb.Fred.home.sl[[col]] <- wb.Fred[[col]] / (sleep_days * 0.5 +  (wb.Fred$time.home / 24 - sleep_days) * sr.fred[col])
}

# 7 days w/sleeping. Pivot wb.Fred.7d.sl to long format
wb.Fred.7d.sl_long <- wb.Fred.7d.sl %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.WB")

# Merge both data frames by the congener column
merged_data <- merge(conc.Fred_long, wb.Fred.7d.sl_long, by = c("ID", "congener"))

# Filter out rows where Conc.Air or Conc.WB are zero
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

filtered_data <- filtered_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
plotfred7dsl <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
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
  scale_fill_manual(values = color_palette) +   # Manual color scale for congener
  scale_shape_manual(values = shape_palette) +  # Manual shape scale for congener
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

plotfred7dsl

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen7dsl.png", plot = plotfred7dsl, width = 5,
       height = 5, dpi = 500)

# Home time. Pivot wb.Fred.7d to long format
wb.Fred.home.sl_long <- wb.Fred.home.sl %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.WB")

# Merge both data frames by the congener column
merged_data <- merge(conc.Fred_long, wb.Fred.home.sl_long, by = c("ID", "congener"))

# Filter out rows where Conc.Air or Conc.WB are zero
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

filtered_data <- filtered_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
plotfredhomesl <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
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
  scale_fill_manual(values = color_palette) +   # Manual color scale for congener
  scale_shape_manual(values = shape_palette) +  # Manual shape scale for congener
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

plotfredhomesl

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksenhomesl.png",
       plot = plotfredhomesl, width = 5, height = 5, dpi = 500)

# Use constant SR & adjust air concentration to exposure time -------------
wb.Fred.7d.v2 <- wb.Fred[, c("ID", "group", matching_cols)]
for (col in matching_cols) {
  wb.Fred.7d.v2[[col]] <- wb.Fred[[col]] / (0.5 * 7) # change SR
}

# 7 days. Pivot wb.Fred.7d to long format
wb.Fred.7d.v2_long <- wb.Fred.7d.v2 %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.WB")

conc.Fred <- data.Fred %>% filter(measurement == "concentration")
# Pivot conc.Fred to long format
conc.Fred_long <- conc.Fred %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.Air")

# Create a copy of conc.Fred to store the modified data
conc.Fred_mod <- conc.Fred
# Identify columns that start with "PCB"
pcb_columns <- grep("^PCB", names(conc.Fred), value = TRUE)
# Adjust the air concentration to the time the volunteer were exposed
# The WB was not 100% of the time exposed to the apartment concentration
conc.Fred_mod[pcb_columns] <- conc.Fred[pcb_columns] * (conc.Fred$time.home / (7 * 24))
# Pivot conc.Fred to long format
conc.Fred_mod_long <- conc.Fred_mod %>%
  pivot_longer(cols = starts_with("PCB"), 
               names_to = "congener", 
               values_to = "Conc.Air")

# Merge both data frames by the congener column
merged_data <- merge(conc.Fred_mod_long, wb.Fred.7d.v2_long, by = c("ID", "congener"))

# Filter out rows where Conc.Air or Conc.WB are zero
filtered_data <- merged_data %>%
  filter(Conc.Air != 0, Conc.WB != 0)

filtered_data <- filtered_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot the data
plotfred7d.v2 <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, 
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
  scale_fill_manual(values = color_palette, labels = names(sr.fred)) +
  scale_shape_manual(values = shape_palette, labels = names(sr.fred)) +
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

plotfred7d.v2

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen7dSrCteTieExpAdj.png", plot = plotfred7d.v2,
       width = 7, height = 7, dpi = 500)


