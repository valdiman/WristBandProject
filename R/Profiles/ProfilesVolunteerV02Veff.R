## Script to fabricate PCB profiles
# and calculate cosine theta for similarity analysis

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("ggpubr")
install.packages("lsa")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(ggpubr) # ggarrange
  library(lsa) # cosine theta function
}

# Read concentration data generated from AirConcVolunteer2Veff.R ----------
{
  conc.air <- read.csv("Output/Data/csv/Volunteer/VolunteerConcStaticWB2.csv")
  conc.wb <- read.csv("Output/Data/csv/Volunteer/VolunteerConcWB2.csv")
}

# Create PCB Profiles -----------------------------------------------------
# (1) Air WBs
tmp.wb.air <- colSums(conc.air[, 2:3], na.rm = TRUE)
prof.air.conc <- sweep(conc.air[, 2:3], 2, tmp.wb.air, FUN = "/")
congener <- conc.air$X
prof.air.conc <- cbind(congener, prof.air.conc)
#Then turn it back into a factor with the levels in the correct order
prof.air.conc$congener <- factor(prof.air.conc$congener,
                                 levels = unique(prof.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.air.conc[, 2:3], na.rm = TRUE)

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
# Vol 1
prof_combined.V1 <- prof.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V1.o) %>%
              rename(Conc = conc.V1.o) %>%
              mutate(Source = "Vol. 1 office-only")) %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V1.h) %>%
              rename(Conc = conc.V1.h) %>%
              mutate(Source = "Vol. 1 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 1 office-only",
                                            "Vol. 1 full-day")))

# Plot
# Create the plot with the legend moved inside
p_prof_comb.V1 <- ggplot(prof_combined.V1, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 1 office-only" = "#009E73",
                               "Vol. 1 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.V1)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Vol1Veff.png",
       plot = p_prof_comb.V1, width = 10, height = 3, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V1 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 1 office-only`, `Vol. 1 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V1 <- ggplot(plot_data, aes(x = `Air PCB office 1`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 1 office-only" = "#009E73",
               "Vol. 1 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V1)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Vol1Veff.png",
       plot = p_scat_comb.V1, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V1 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 1` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V1

# Vol 2
prof_combined.V2 <- prof.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V2.o) %>%
              rename(Conc = conc.V2.o) %>%
              mutate(Source = "Vol. 2 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V2.h) %>%
              rename(Conc = conc.V2.h) %>%
              mutate(Source = "Vol. 2 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 2 office-only",
                                            "Vol. 2 full-day")))

# Plots
p_prof_comb.V2 <- ggplot(prof_combined.V2, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 2 office-only" = "#009E73",
                               "Vol. 2 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.V2)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Vol2Veff.png",
       plot = p_prof_comb.V2, width = 10, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V2 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 2 office-only`, `Vol. 2 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V2 <- ggplot(plot_data, aes(x = `Air PCB office 1`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 2 office-only" = "#009E73",
               "Vol. 2 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V2)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Vol2Veff.png",
       plot = p_scat_comb.V2, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V2 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 1` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V2

# Vol 3
prof_combined.V3 <- prof.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V3.o) %>%
              rename(Conc = conc.V3.o) %>%
              mutate(Source = "Vol. 3 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V3.h) %>%
              rename(Conc = conc.V3.h) %>%
              mutate(Source = "Vol. 3 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 3 office-only",
                                            "Vol. 3 full-day")))

# Plots
p_prof_comb.V3 <- ggplot(prof_combined.V3, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 3 office-only" = "#009E73",
                               "Vol. 3 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.V3)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Vol3Veff.png",
       plot = p_prof_comb.V3, width = 10, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V3 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 3 office-only`, `Vol. 3 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V3 <- ggplot(plot_data, aes(x = `Air PCB office 1`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 3 office-only" = "#009E73",
               "Vol. 3 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V3)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Vol3Veff.png",
       plot = p_scat_comb.V3, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V3 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 1` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V3

# Vol 8 (data say Vol 5)
prof_combined.V8 <- prof.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V5.o) %>%
              rename(Conc = conc.V5.o) %>%
              mutate(Source = "Vol. 8 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V5.h) %>%
              rename(Conc = conc.V5.h) %>%
              mutate(Source = "Vol. 8 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 8 office-only",
                                            "Vol. 8 full-day")))

# Plots
p_prof_comb.V8 <-  ggplot(prof_combined.V8, aes(x = congener, y = Conc,
                                                fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 8 office-only" = "#009E73",
                               "Vol. 8 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.V8)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Vol8Veff.png",
       plot = p_prof_comb.V8, width = 10, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V8 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 8 office-only`, `Vol. 8 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V8 <- ggplot(plot_data, aes(x = `Air PCB office 2`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 8 office-only" = "#009E73",
               "Vol. 8 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V8)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Vol8Veff.png",
       plot = p_scat_comb.V8, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V8 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 2` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V8

# Vol 9 (data say Vol 4)
prof_combined.V9 <- prof.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V4.o) %>%
              rename(Conc = conc.V4.o) %>%
              mutate(Source = "Vol. 9 office-only")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V4.h) %>%
              rename(Conc = conc.V4.h) %>%
              mutate(Source = "Vol. 9 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 9 office-only",
                                            "Vol. 9 full-day")))

# Plots
p_prof_comb.V9 <- ggplot(prof_combined.V9, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 9 office-only" = "#009E73",
                               "Vol. 9 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.V9)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Vol9Veff.png",
       plot = p_prof_comb.V9, width = 10, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V9 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 9 office-only`, `Vol. 9 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V9 <- ggplot(plot_data, aes(x = `Air PCB office 2`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 9 office-only" = "#009E73",
               "Vol. 9 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V9)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Vol9Veff.png",
       plot = p_scat_comb.V9, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V9 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 2` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V9

# Volunteers in location 1 (Conc.Air.1) Vol. 2 out
prof_combined.1 <-  prof.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V3.h) %>%
              rename(Conc = conc.V3.h) %>%
              mutate(Source = "Vol. 3 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V1.h) %>%
              rename(Conc = conc.V1.h) %>%
              mutate(Source = "Vol. 1 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 1 full-day",
                                            "Vol. 3 full-day")))

p_prof_comb.1 <- ggplot(prof_combined.1, aes(x = congener, y = Conc,
                                              fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 1 full-day" = "#009E73",
                               "Vol. 3 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.1)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Office1Veff.png",
       plot = p_prof_comb.1, width = 10, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.1 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 1 full-day`, `Vol. 3 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.1 <- ggplot(plot_data, aes(x = `Air PCB office 1`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 1 full-day" = "#009E73",
               "Vol. 3 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.1)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Office1Veff.png",
       plot = p_scat_comb.1, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.1 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 1` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.1

# Volunteers in location 2 (Conc.Air.2)
prof_combined.2 <-  prof.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V4.h) %>%
              rename(Conc = conc.V4.h) %>%
              mutate(Source = "Vol. 8 full-day")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V5.h) %>%
              rename(Conc = conc.V5.h) %>%
              mutate(Source = "Vol. 9 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 8 full-day",
                                            "Vol. 9 full-day")))

p_prof_comb.2 <- ggplot(prof_combined.2, aes(x = congener, y = Conc,
                                             fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
           linewidth = 0.2) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  ylab(expression(bold("Conc. Fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 8 full-day" = "#009E73",
                               "Vol. 9 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))

# Print the plots
print(p_prof_comb.2)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Barplot/prof_combined.Office2Veff.png",
       plot = p_prof_comb.2, width = 10, height = 3, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.2 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 8 full-day`, `Vol. 9 full-day`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.2 <- ggplot(plot_data, aes(x = `Air PCB office 2`, y = Conc, color = Vol_Type)) +
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
    values = c("Vol. 8 full-day" = "#009E73",
               "Vol. 9 full-day" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.2)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/Scatterplot/prof_combined.Office2Veff.png",
       plot = p_scat_comb.2, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.2 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office 2` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.2

# Calculate cosine theta --------------------------------------------------
# Need to change the format of prof_combined...
# Create a term-document matrix from the 'Source' variable
# Vol 1
prof_combined_wide.V1  <- prof_combined.V1 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V1[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V1 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V1

# Vol 2
prof_combined_wide.V2  <- prof_combined.V2 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V2[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V2 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V2

# Vol 3
prof_combined_wide.V3  <- prof_combined.V3 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V3[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V3 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V3

# Vol 8
prof_combined_wide.V8  <- prof_combined.V8 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V8[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V8 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V8

# Vol 9
prof_combined_wide.V9  <- prof_combined.V9 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V9[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V9 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V9

# Between volunteers from the same office
# Room/office 1
prof_combined_wide.1  <- prof_combined.1 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.1[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.1 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.1

# Room/office 2
prof_combined_wide.2  <- prof_combined.2 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.2[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.2 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.2



