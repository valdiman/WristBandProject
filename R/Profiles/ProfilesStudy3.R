## Script to fabricate PCB profiles
# and calculate cosine theta for similarity analysis
# for Study 3

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

# Read concentration data generated from AirConcVolunteerVeffStudy3.R -----------
{
  conc.air <- read.csv("Output/Data/Study3/VolunteerConcStaticWBStudy3.csv")
  conc.wb <- read.csv("Output/Data/Study3/VolunteerConcWBStudy3.csv")
}

# Have the same number of congeners
common_ids <- intersect(conc.air$X, conc.wb$X)
conc.air <- conc.air[conc.air$X %in% common_ids, ]
conc.wb  <- conc.wb[conc.wb$X %in% common_ids, ]

# Create PCB profiles -----------------------------------------------------
# (1) Air WBs
tmp.wb.air <- colSums(conc.air[, 2:8], na.rm = TRUE)
prof.air.conc <- sweep(conc.air[, 2:8], 2, tmp.wb.air, FUN = "/")
congener <- conc.air$X
prof.air.conc <- cbind(congener, prof.air.conc)
#Then turn it back into a factor with the levels in the correct order
prof.air.conc$congener <- factor(prof.air.conc$congener,
                                levels = unique(prof.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.air.conc[, 2:8], na.rm = TRUE)

# (2) Wore WBs
tmp.wb.wr <- colSums(conc.wb[, 2:16], na.rm = TRUE)
prof.wb.conc <- sweep(conc.wb[, 2:16], 2, tmp.wb.wr, FUN = "/")
prof.wb.conc <- cbind(congener, prof.wb.conc)
# Then turn it back into a factor with the levels in the correct order
prof.wb.conc$congener <- factor(prof.wb.conc$congener,
                               levels = unique(prof.wb.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.conc[, 2:16], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# V4
prof_combined.V4 <- prof.air.conc %>%
  select(congener, Conc.Air.V4) %>%
  rename(Conc = Conc.Air.V4) %>%
  mutate(Source = "Air PCB office") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V4.nd) %>%
              rename(Conc = Conc.WB.V4.nd) %>%
              mutate(Source = "V4 nd")) %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V4.d) %>%
              rename(Conc = Conc.WB.V4.d) %>%
              mutate(Source = "V4 d"))  # Change this label

# Zoom plot (Figure in manuscript)
prof_combined.V4.z <- prof_combined.V4 %>%
  mutate(
    congener_chr = as.character(congener),
    congener_num = as.numeric(str_extract(congener_chr, "\\d+\\.?\\d*"))
  )

prof_subset <- prof_combined.V4.z %>%
  filter(congener_num >= 42 & congener_num <= 120) %>%
  mutate(congener_chr = gsub("\\.", "+", congener_chr))

p_prof_comb.V4_subset <- ggplot(prof_subset, aes(x = congener, y = Conc, fill = Source)) +
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
                               "V4 nd" = "#009E73",
                               "V4 d" = "#E69F00"),
                    guide = guide_legend()) +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# See plot
p_prof_comb.V4_subset

ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V4VeffZoom.png",
       plot = p_prof_comb.V4_subset, width = 10, height = 3, dpi = 500)

# Bar plot all congeners
# Create 3 plots
p_prof_comb.V4 <- ggplot(prof_combined.V4, aes(x = congener,
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
    legend.position = "none") + #,
    # To remove congeners in x-axis
    #axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
    #                                  size = 9, face = "bold"),
    #axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue", "V4 nd" = "#009E73",
                               "V4 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V4)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V4VeffWOCongeners.png",
       plot = p_prof_comb.V4, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V4 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V4 nd`, `V4 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V4 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V4 nd" = "#009E73",
               "V4 d" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V4)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V4Veff.png",
       plot = p_scat_comb.V4, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V4 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V4

# V5
prof_combined.V5 <- prof.air.conc %>%
  select(congener, Conc.Air.V5) %>%
  rename(Conc = Conc.Air.V5) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V5.d) %>%
              rename(Conc = Conc.WB.V5.d) %>%
              mutate(Source = "V5 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V5.nd) %>%
              rename(Conc = Conc.WB.V5.nd) %>%
              mutate(Source = "V5 nd"))

# Plots
# Bar plot all congeners
# Create 3 plots
p_prof_comb.V5 <- ggplot(prof_combined.V5, aes(x = congener,
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
    legend.position = "none") + #,
  # To remove congeners in x-axis
  #axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1,
  #                                  size = 9, face = "bold"),
  #axis.ticks.x.bottom = element_line()) +
  scale_fill_manual(values = c("Air PCB office" = "blue", "V5 nd" = "#009E73",
                               "V5 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See the plots
print(p_prof_comb.V5)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V5VeffWOCongeners.png",
       plot = p_prof_comb.V5, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V5 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V5 nd`, `V5 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V5 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V5 nd" = "#009E73",
               "V5 d" = "#E69F00"),
    guide = guide_legend())

# See the plots
print(p_scat_comb.V5)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V5Veff.png",
       plot = p_scat_comb.V5, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V5 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V5

# V6
prof_combined.V6 <- prof.air.conc %>%
  select(congener, Conc.Air.V6) %>%
  rename(Conc = Conc.Air.V6) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V6.nd) %>%
              rename(Conc = Conc.WB.V6.nd) %>%
              mutate(Source = "V6 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V6.d) %>%
              rename(Conc = Conc.WB.V6.d) %>%
              mutate(Source = "V6 d"))

# Create 3 plots
p_prof_comb.V6 <- ggplot(prof_combined.V6, aes(x = congener,
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
  scale_fill_manual(values = c("Air PCB office" = "blue", "V6 nd" = "#009E73",
                               "V6 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See the plots
print(p_prof_comb.V6)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V6Veff.png",
       plot = p_prof_comb.V6, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V6 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V6 nd`, `V6 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V6 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V6 nd" = "#009E73",
               "V6 d" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V6)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V6Veff.png",
       plot = p_scat_comb.V6, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V6 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V6

# V7
prof_combined.V7 <- prof.air.conc %>%
  select(congener, Conc.Air.V7) %>%
  rename(Conc = Conc.Air.V7) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V7.nd) %>%
              rename(Conc = Conc.WB.V7.nd) %>%
              mutate(Source = "V7 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V7.d) %>%
              rename(Conc = Conc.WB.V7.d) %>%
              mutate(Source = "V7 d"))

# Create 3 plots
p_prof_comb.V7 <- ggplot(prof_combined.V7, aes(x = congener,
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
  scale_fill_manual(values = c("Air PCB office" = "blue", "V7 nd" = "#009E73",
                               "V7 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V7)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V7Veff.png",
       plot = p_prof_comb.V7, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V7 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V7 nd`, `V7 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V7 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V7 nd" = "#009E73",
               "V7 d" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V7)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V7Veff.png",
       plot = p_scat_comb.V7, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V7 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V7

# V8
prof_combined.V8 <- prof.air.conc %>%
  select(congener, Conc.Air.V8) %>%
  rename(Conc = Conc.Air.V8) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V8.d) %>%
              rename(Conc = Conc.WB.V8.d) %>%
              mutate(Source = "V8 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V8.nd) %>%
              rename(Conc = Conc.WB.V8.nd) %>%
              mutate(Source = "V8 nd"))

# Create 3 plots
p_prof_comb.V8 <- ggplot(prof_combined.V8, aes(x = congener,
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
  scale_fill_manual(values = c("Air PCB office" = "blue", "V8 nd" = "#009E73",
                               "V8 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V8)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V8Veff.png",
       plot = p_prof_comb.V8, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V8 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V8 nd`, `V8 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V8 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V8 nd" = "#009E73",
               "V8 d" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V8)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V8Veff.png",
       plot = p_scat_comb.V8, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V8 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V8

# V9
prof_combined.V9 <- prof.air.conc %>%
  select(congener, Conc.Air.V9) %>%
  rename(Conc = Conc.Air.V9) %>%
  mutate(Source = "Air PCB office") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V9.nd) %>%
              rename(Conc = Conc.WB.V9.nd) %>%
              mutate(Source = "V9 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, Conc.WB.V9.d) %>%
              rename(Conc = Conc.WB.V9.d) %>%
              mutate(Source = "V9 d"))

# Create 3 plots
p_prof_comb.V9 <- ggplot(prof_combined.V9, aes(x = congener,
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
  scale_fill_manual(values = c("Air PCB office" = "blue", "V9 nd" = "#009E73",
                               "V9 d" = "#E69F00")) + 
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# Print the plots
print(p_prof_comb.V9)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Barplot/prof_combined.V9Veff.png",
       plot = p_prof_comb.V9, width = 22, height = 5, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V9 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`V9 nd`, `V9 d`),
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V9 <- ggplot(plot_data, aes(x = `Air PCB office`, y = Conc, color = Vol_Type)) +
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
    legend.position = c(0.15, 0.9)) +
  scale_color_manual(
    values = c("V9 nd" = "#009E73",
               "V9 d" = "#E69F00"),
    guide = guide_legend())

# Print the plots
print(p_scat_comb.V9)

# Save plot in folder
ggsave("Output/Plots/Profiles/Study3/Scatterplot/prof_combined.V9Veff.png",
       plot = p_scat_comb.V9, width = 5, height = 5, dpi = 500)

# See max abs difference
diff.V9 <- plot_data %>%
  mutate(abs_diff = abs(`Air PCB office` - Conc)) %>%
  slice_max(abs_diff, n = 1, with_ties = FALSE)
diff.V9

# Calculate cosine theta --------------------------------------------------
# Need to change the format of prof_combined...
# Create a term-document matrix from the 'Source' variable
# V4
prof_combined_wide.V4  <- prof_combined.V4 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V4[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V4 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V4

# V5
prof_combined_wide.V5  <- prof_combined.V5 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V5[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V5 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V5

# V6
prof_combined_wide.V6  <- prof_combined.V6 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V6[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V6 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V6

# V7
prof_combined_wide.V7  <- prof_combined.V7 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V7[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V7 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V7

# V8
prof_combined_wide.V8  <- prof_combined.V8 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V8[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V8 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V8

# V9
prof_combined_wide.V9  <- prof_combined.V9 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V9[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V9 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V9

