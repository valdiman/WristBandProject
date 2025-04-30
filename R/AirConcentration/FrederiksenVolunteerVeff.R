## Script to combine Frederiksen and Volunteer 2 data to be plotted

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
  vol2 <- read.csv("Output/Data/csv/FrederiksenPCB/Volunteer2_PCBiVeff.csv")
  fred <- read.csv("Output/Data/csv/FrederiksenPCB/Frederiksen_PCBiVeff.csv")
}

# Format data -------------------------------------------------------------
fred_trimmed <- fred %>%
  rename(Conc.WB = est_conc,
         Conc.Air = obs_conc) %>%
  select(congener, Conc.Air, Conc.WB, ID) %>%
  mutate(source = "fred")

vol2 <- vol2 %>%
  mutate(source = "vol2")

combined_data <- bind_rows(fred_trimmed, vol2)

# Create a custom order for congener_names
congener_names <- c("PCB8", "PCB18.30", "PCB20.28", "PCB31", "PCB44.47.65",
                    "PCB52", "PCB66", "PCB99", "PCB90.101.113", "PCB105",
                    "PCB118", "PCB129.138.163", "PCB153.168", "PCB180.193")

# Apply the factor levels to the congener column in the dataset
combined_data <- combined_data %>%
  mutate(congener = factor(congener, levels = congener_names))

combined_data$source <- factor(combined_data$source)

# Define the custom color palette
color_palette2 <- c("red", "blue", "green", "yellow", "purple", "orange", 
                    "pink", "cyan", "brown", "magenta", "limegreen", 
                    "skyblue", "gold", "slategray")  # 14 colors as an example

shape_palette2 <- c(21, 24)  # Shape 21 for one source and 24 for another

# Plot
plot.FerdVol2 <- ggplot(combined_data, aes(x = Conc.Air, y = Conc.WB, fill = congener, shape = source)) +
  geom_point(size = 2, stroke = 0.5, color = "black") +  # Black borders, solid color fill
  theme_bw() +
  theme(aspect.ratio = 1) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(0.001, 10^4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^4),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Indoor Air Concentration PCBi (ng/m"^3*")"))) +
  ylab(expression(bold("Predicted Concentration PCBi (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  scale_shape_manual(values = shape_palette2, labels = c("Frederiksen et. al. 2022",
                                                         "This study")) +
  scale_fill_manual(values = color_palette2) +  # Apply custom color palette for congener
  labs(fill = "Congener", shape = "Data Source") +
  theme(legend.position = "right") + 
  guides(fill = guide_legend(override.aes = list(shape = 21)))

# See plot
plot.FerdVol2

# Save plot in folder
ggsave("Output/Plots/AirConcentrations/Frederiksen/FrederiksenVol2PCBiVeff.png",
       plot = plot.FerdVol2, width = 6, height = 6, dpi = 500)

