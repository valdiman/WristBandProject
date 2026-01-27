## Script to compare mass accumulated in WBs and air concentration
# from Frederiksen et al. (2022).
# Airborne concentrations of individual PCB congeners or co-eluting congeners
# measured in apartments (active low volume sampler), and masses of individual
# PCB congeners or co-eluting congeners accumulated in wristbands worn by
# participants

# Install packages
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("RColorBrewer")

# Load libraries
{
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
  library(scales)
}

# Read data ---------------------------------------------------------------
{
  data.Fred.mass <- read.csv("Data/IRO/10_FrederiksenWBMass2022.csv",
                             check.names = FALSE)
  data.Fred.conc <- read.csv("Data/IRO/11_FrederiksenConc2022.csv",
                             check.names = FALSE)
  # Read individual PCB logKoa
  logKoa <- read.csv("Data/IRO/12_logKoa.csv")
  # Read ko generated from PersonalSamplingRates.R file
  ko.p <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSR.csv")
}

# Select data -------------------------------------------------------------
congener_names <- grep("^PCB", colnames(data.Fred.mass), value = TRUE)
# Select ko values
ko.p <- ko.p %>% select(congener, Average_ko)
# Filter ko to keep only matching congeners
ko_fred <- ko.p[ko.p$congener %in% congener_names, ]
# Extract only the ko values
ko.fred <- ko_fred$Average_ko
# Name the vector with congener names for reference
names(ko.fred) <- ko_fred$congener

# Calculate logKwb --------------------------------------------------------
# Regression created with data from Tromp et al 2019 (Table 2, Wristband)
# & Frederiksen et al 2022 (Table 3)
{
  logKwb <- data.frame(
    congener = logKoa$congener,
    logKwb = 0.6156 * logKoa$logKoa + 2.161) # R2 = 0.96
  logKwb_fred <- logKwb$logKwb[logKwb$congener %in% congener_names]
}

# Calculate Veffs only with home time -------------------------------------
# Change time to days
wb.mass.Fred <- data.Fred.mass %>%
  mutate(time = time / 24)

# Define Vwb and Awb
Vwb <- 0.00000473 # [m3]
Awb <- 0.0054773 # [m2]

# Veff calculations
veff_fred <- outer(
  wb.mass.Fred$time,
  1:14,
  Vectorize(function(t, i) {
    10^(logKwb_fred)[i] * Vwb * (1 - exp(-ko_fred$Average_ko[i] * Awb / Vwb / 10^(logKwb_fred)[i] * t))
  })
)

# Add row and column names
colnames(veff_fred) <- ko_fred$congener
rownames(veff_fred) <- wb.mass.Fred$sid

# Extract PCB columns by name from wb.mass.Fred using column names of veff_fred
pcb_mass <- wb.mass.Fred[, congener_names]

# Set row names to match sid
rownames(pcb_mass) <- wb.mass.Fred$sid

# Ensure column order matches veff_fred
pcb_mass <- pcb_mass[, colnames(veff_fred)]

# Calculate concentration
conc_fred <- pcb_mass / veff_fred

conc_fred_long <- as.data.frame(conc_fred) %>%
  mutate(sid = rownames(conc_fred)) %>%
  pivot_longer(-sid, names_to = "congener", values_to = "est_conc")

# Filter data.Fred.conc to keep only "concentration" rows
data_conc_long <- data.Fred.conc %>%
  select(sid, starts_with("PCB")) %>%
  pivot_longer(-sid, names_to = "congener", values_to = "obs_conc")

# Join the datasets by ID and congener
compare_df <- left_join(conc_fred_long, data_conc_long, by = c("sid", "congener"))

plot_data <- compare_df %>%
  filter(est_conc > 0, obs_conc > 0)

plot_data <- plot_data %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot
plotfred <- ggplot(plot_data, aes(x = obs_conc, y = est_conc, 
                                      fill = congener, shape = congener)) +
  geom_point(size = 2.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(0.001, 1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Apartment Air Concentration PCBi (ng/m"^3*")"))) +
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
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

# Show the plot
plotfred

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen/FrederiksenHomeVeff.png", plot = plotfred,
       width = 7, height = 7, dpi = 500)

# Export data
write.csv(plot_data,
          file = "Output/Data/csv/FrederiksenPCB/Frederiksen_PCBiVeff.csv")

# Calculate Veff for 7 days -----------------------------------------------
# Change time to 7 days
wb.mass.Fred.2 <- data.Fred.mass %>%
  mutate(time = 7 * 24)

# Veff calculations
veff_fred.2 <- outer(
  wb.mass.Fred.2$time,
  1:14,
  Vectorize(function(t, i) {
    10^(logKwb_fred)[i] * Vwb * (1 - exp(-ko_fred$Average_ko[i] * Awb / Vwb / 10^(logKwb_fred)[i] * t))
  })
)

# Add row and column names
colnames(veff_fred.2) <- ko_fred$congener
rownames(veff_fred.2) <- wb.mass.Fred.2$sid

# Extract PCB columns by name from wb.mass.Fred using column names of veff_fred
pcb_mass <- wb.mass.Fred.2[, congener_names]

# Set row names to match IDs
rownames(pcb_mass) <- wb.mass.Fred.2$sid

# Ensure column order matches veff_fred.2
pcb_mass <- pcb_mass[, colnames(veff_fred.2)]

# Element-wise division
conc_fred.2 <- pcb_mass / veff_fred.2

conc_fred_long.2 <- as.data.frame(conc_fred.2) %>%
  mutate(sid = rownames(conc_fred.2)) %>%
  pivot_longer(-sid, names_to = "congener", values_to = "est_conc")

# Filter data.Fred.conc to keep only "concentration" rows
data_conc_long <- data.Fred.conc %>%
  select(sid, starts_with("PCB")) %>%
  pivot_longer(-sid, names_to = "congener", values_to = "obs_conc")

# Join the datasets by ID and congener
compare_df.2 <- left_join(conc_fred_long.2, data_conc_long, by = c("sid", "congener"))

plot_data.2 <- compare_df.2 %>%
  filter(est_conc > 0, obs_conc > 0)

plot_data.2 <- plot_data.2 %>%
  mutate(congener = factor(congener, levels = congener_names))

color_palette <- c("red", "blue", "green", "purple", "orange", "brown", 
                   "pink", "yellow", "cyan", "gray", "black", "violet", 
                   "magenta", "indianred") # 14 colors as an example

shape_palette <- c(21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 21, 21, 22, 22)

# Plot
plotfred.2 <- ggplot(plot_data.2, aes(x = obs_conc, y = est_conc, 
                                  fill = congener, shape = congener)) +
  geom_point(size = 2.5, color = "black", stroke = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(0.0001, 1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.0001, 1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Apartment Air Concentration PCBi (ng/m"^3*")"))) +
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
  labs(fill = "Congener", shape = "Congener") +
  theme(legend.position = "right")

# Show the plot
plotfred.2

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen/Frederiksen7dVeff.png", plot = plotfred.2,
       width = 7, height = 7, dpi = 500)

# Export data
write.csv(plot_data.2,
          file = "Output/Data/csv/FrederiksenPCB/Frederiksen_PCBiV3.csv")
