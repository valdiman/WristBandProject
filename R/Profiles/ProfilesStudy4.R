## Script to fabricate PCB profiles
# and calculate cosine theta for similarity analysis
# for Study 4

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("stringr")
install.packages("ggpubr")
install.packages("lsa")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggpubr) # ggarrange
  library(lsa) # cosine theta function
}

# Read concentration data generated from AirConcVolunteerVeffStudy4.R ----------
{
  conc.air <- read.csv("Output/Data/Study4/VolunteerConcStaticWBStudy4.csv")
  conc.wb <- read.csv("Output/Data/Study4/VolunteerConcWBStudy4.csv")
}

# Create PCB Profiles -----------------------------------------------------
# (1) Air WBs
tmp.wb.air <- colSums(conc.air[, 2:6], na.rm = TRUE)
prof.air.conc <- sweep(conc.air[, 2:6], 2, tmp.wb.air, FUN = "/")
congener <- conc.air$X
prof.air.conc <- cbind(congener, prof.air.conc)
#Then turn it back into a factor with the levels in the correct order
prof.air.conc$congener <- factor(prof.air.conc$congener,
                                 levels = unique(prof.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.air.conc[, 2:6], na.rm = TRUE)

# (2) Wore WBs
tmp.wb.wr <- colSums(conc.wb[, 2:11], na.rm = TRUE)
prof.wb.conc <- sweep(conc.wb[, 2:11], 2, tmp.wb.wr, FUN = "/")
prof.wb.conc <- cbind(congener, prof.wb.conc)
#Then turn it back into a factor with the levels in the correct order
prof.wb.conc$congener <- factor(prof.wb.conc$congener,
                                levels = unique(prof.wb.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.conc[, 2:11], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# Vol 15
prof_combined.V15 <- prof.air.conc %>%
  select(congener, Conc.Air.V15) %>%
  rename(Conc = Conc.Air.V15) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V15.o) %>%
              rename(Conc = Conc.WB.V15.o) %>%
              mutate(Source = "V15 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V15.d) %>%
              rename(Conc = Conc.WB.V15.d) %>%
              mutate(Source = "V15 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V15 office-only",
                                            "V15 full-day")))

# Create 3 plots
p_prof_comb.V15 <- ggplot(prof_combined.V15, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V15 office-only" = "#009E73",
                               "V15 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V15)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.V15Veff.png",
       plot = p_prof_comb.V15, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V15 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V15 office-only`, `V15 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V15 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V15 office-only" = "#009E73",
               "V15 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V15)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.Vol15Veff.png",
       plot = p_scat_comb.V15, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V15 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V15

# Vol 16
prof_combined.V16 <- prof.air.conc %>%
  select(congener, Conc.Air.V16) %>%
  rename(Conc = Conc.Air.V16) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V16.o) %>%
              rename(Conc = Conc.WB.V16.o) %>%
              mutate(Source = "V16 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V16.d) %>%
              rename(Conc = Conc.WB.V16.d) %>%
              mutate(Source = "V16 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V16 office-only",
                                            "V16 full-day")))

# Create 3 plots
p_prof_comb.V16 <- ggplot(prof_combined.V16, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V16 office-only" = "#009E73",
                               "V16 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V16)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.V16Veff.png",
       plot = p_prof_comb.V16, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V16 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V16 office-only`, `V16 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V16 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V16 office-only" = "#009E73",
               "V16 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V16)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.Vol16Veff.png",
       plot = p_scat_comb.V16, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V16 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V16

# Vol 17
prof_combined.V17 <- prof.air.conc %>%
  select(congener, Conc.Air.V17) %>%
  rename(Conc = Conc.Air.V17) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V17.o) %>%
              rename(Conc = Conc.WB.V17.o) %>%
              mutate(Source = "V17 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V17.d) %>%
              rename(Conc = Conc.WB.V17.d) %>%
              mutate(Source = "V17 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V17 office-only",
                                            "V17 full-day")))

# Create 3 plots
p_prof_comb.V17 <- ggplot(prof_combined.V17, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V17 office-only" = "#009E73",
                               "V17 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V17)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.V17Veff.png",
       plot = p_prof_comb.V17, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V17 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V17 office-only`, `V17 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V17 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V17 office-only" = "#009E73",
               "V17 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V17)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.Vol17Veff.png",
       plot = p_scat_comb.V17, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V17 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V17

# Vol 18
prof_combined.V18 <- prof.air.conc %>%
  select(congener, Conc.Air.V18) %>%
  rename(Conc = Conc.Air.V18) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V18.o) %>%
              rename(Conc = Conc.WB.V18.o) %>%
              mutate(Source = "V18 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V18.d) %>%
              rename(Conc = Conc.WB.V18.d) %>%
              mutate(Source = "V18 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V18 office-only",
                                            "V18 full-day")))

# Create 3 plots
p_prof_comb.V18 <- ggplot(prof_combined.V18, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V18 office-only" = "#009E73",
                               "V18 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V18)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.Vol18Veff.png",
       plot = p_prof_comb.V18, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V18 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V18 office-only`, `V18 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V18 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V18 office-only" = "#009E73",
               "V18 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V18)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.Vol18Veff.png",
       plot = p_scat_comb.V18, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V18 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V18

# Vol 19
prof_combined.V19 <- prof.air.conc %>%
  select(congener, Conc.Air.V19) %>%
  rename(Conc = Conc.Air.V19) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V19.o) %>%
              rename(Conc = Conc.WB.V19.o) %>%
              mutate(Source = "V19 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V19.d) %>%
              rename(Conc = Conc.WB.V19.d) %>%
              mutate(Source = "V19 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V19 office-only",
                                            "V19 full-day")))

# Create 3 plots
p_prof_comb.V19 <- ggplot(prof_combined.V19, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V19 office-only" = "#009E73",
                               "V19 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V19)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.Vol19Veff.png",
       plot = p_prof_comb.V19, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V19 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V19 office-only`, `V19 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V19 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V19 office-only" = "#009E73",
               "V19 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V19)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.Vol19Veff.png",
       plot = p_scat_comb.V19, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V19 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V19

# Volunteers V15 and V17 in same office
prof_combined.V15_17 <-  prof.air.conc %>%
  select(congener, Conc.Air.V15) %>%
  rename(Conc = Conc.Air.V15) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V17.d) %>%
              rename(Conc = Conc.WB.V17.d) %>%
              mutate(Source = "V17 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V15.d) %>%
              rename(Conc = Conc.WB.V15.d) %>%
              mutate(Source = "V15 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V15 full-day",
                                            "V17 full-day")))

# Create 3 plots
p_prof_comb.15_17 <- ggplot(prof_combined.V15_17, aes(x = congener,
                                               y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V15 full-day" = "#009E73",
                               "V17 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.15_17)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.fulldayV15_17Veff.png",
       plot = p_prof_comb.15_17, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V15_17 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V15 full-day`, `V17 full-day`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.15_17 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V15 full-day" = "#009E73",
               "V17 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.15_17)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.V15_17Veff.png",
       plot = p_scat_comb.15_17, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.15_17 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.15_17

# Volunteers 18 and 19, same office
prof_combined.18_19 <-  prof.air.conc %>%
  select(congener, Conc.Air.V18) %>%
  rename(Conc = Conc.Air.V18) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V18.d) %>%
              rename(Conc = Conc.WB.V18.d) %>%
              mutate(Source = "V18 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V19.d) %>%
              rename(Conc = Conc.WB.V19.d) %>%
              mutate(Source = "V19 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V18 full-day",
                                            "V19 full-day")))

# Zoom plot (Figure in manuscript)
prof_combined.18_19.z <- prof_combined.18_19 %>%
  mutate(
    congener_chr = as.character(congener),
    congener_num = as.numeric(str_extract(congener_chr, "\\d+\\.?\\d*"))
  )

prof_subset <- prof_combined.18_19.z %>%
  filter(congener_num >= 42 & congener_num <= 120) %>%
  mutate(congener_chr = gsub("\\.", "+", congener_chr))

p_prof_comb.2_subset <- ggplot(prof_subset, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black", linewidth = 0.2) +
  scale_x_discrete(
    labels = prof_subset$congener_chr) +
  xlab("") + ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20,
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          size = 8,
          face = "bold"),
        axis.ticks.x = element_line(),
        legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V18 full-day" = "#009E73",
                               "V19 full-day" = "#E69F00"),
                    guide = guide_legend()) +
  annotate("text", x = -Inf, y = Inf, label = "(c)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

p_prof_comb.2_subset

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.OfficeV18_19VeffZoom.png",
       plot = p_prof_comb.2_subset, width = 10, height = 5, dpi = 500)

# Create 3 plots
p_prof_comb.2 <- ggplot(prof_combined.2, aes(x = congener,
                                             y = Conc, fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous( limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Conc. Fraction " *Sigma*"PCB"))) +
  theme(
    aspect.ratio = NULL,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                      size = 9, face = "bold"),
    axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue",
                               "V18 full-day" = "#009E73",
                               "V19 full-day" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.2)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Barplot/prof_combined.fulldayV18_19Veff.png",
       plot = p_prof_comb.2, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.2 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V18 full-day`, `V19 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.2 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  xlab(expression(bold("Air Conc. Fraction " *Sigma*"PCB"))) +
  ylab(expression(bold("Volunteer Predited Air Conc. Fraction " *Sigma*"PCB"))) +
  ylim(0, 0.15) +
  xlim(0, 0.15) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    legend.position = c(0.2, 0.9)) +
  scale_color_manual(
    values = c("V18 full-day" = "#009E73",
               "V19 full-day" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.2)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study4/Scatterplot/prof_combined.OfficeV18_19Veff.png",
       plot = p_scat_comb.2, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.2 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.2

# Calculate cosine theta --------------------------------------------------
# Need to change the format of prof_combined...
# Create a term-document matrix from the 'Source' variable
# Vol 15
prof_combined_wide.V15  <- prof_combined.V15 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V15[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V15 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V15

# Vol 16
prof_combined_wide.V16  <- prof_combined.V16 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V16[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V16 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V16

# Vol 17
prof_combined_wide.V17  <- prof_combined.V17 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V17[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V17 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V17

# Vol 18
prof_combined_wide.V18  <- prof_combined.V18 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V18[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V18 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V18

# Vol 19
prof_combined_wide.V19  <- prof_combined.V19 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V19[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V19 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V19

# Between volunteers from the same office (V15, 16, 17)
# Room/office 1
prof_combined.V15_16_17 <- prof.air.conc %>%
  select(congener, Conc.Air.V15) %>%
  rename(Conc = Conc.Air.V15) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V15.d) %>%
              rename(Conc = Conc.WB.V15.d) %>%
              mutate(Source = "V15 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V16.d) %>%
              rename(Conc = Conc.WB.V16.d) %>%
              mutate(Source = "V16 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V17.d) %>%
              rename(Conc = Conc.WB.V17.d) %>%
              mutate(Source = "V17 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office", "V15 full-day",
                                            "V16 full-day", "V17 full-day")))

prof_combined_wide.1  <- prof_combined.V15_16_17 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.1[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Office1 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Office1

# Room/office 2
prof_combined_wide.2  <- prof_combined.18_19 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.2[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.2 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.2

