## Script to fabricate PCB profiles &
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

# Read concentration data generated from AirConcVolunteerVeff.R -----------
{
  conc.air <- read.csv("Output/Data/csv/Volunteer/VolunteerConcStaticWB.csv")
  conc.wb <- read.csv("Output/Data/csv/Volunteer/VolunteerConcWB.csv")
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
# Vol 1
prof_combined.V1 <- prof.air.conc %>%
  select(congener, Conc.Air.V1.V2) %>%
  rename(Conc = Conc.Air.V1.V2) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V1.l) %>%
              rename(Conc = conc.V1.l) %>%
              mutate(Source = "Vol. 1 nd")) %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V1.r) %>%
              rename(Conc = conc.V1.r) %>%
              mutate(Source = "Vol. 1 d"))  # Change this label

# Bar plot
# Create the plot with the legend moved inside
p_prof_comb.V1 <- ggplot(prof_combined.V1, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 1 nd" = "#009E73",
                               "Vol. 1 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(a)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V1)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol1VeffV2.png",
       plot = p_prof_comb.V1, width = 10, height = 3, dpi = 500)

# Scatter 1:1 plot
# Convert to wide format
prof_wide <- prof_combined.V1 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 1 nd`, `Vol. 1 d`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

p_scat_comb.V1 <- ggplot(plot_data, aes(x = `Air PCB`, y = Conc, color = Vol_Type)) +
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
    axis.title.x = element_text(face = "bold", size = 13)) +
  scale_color_manual(
    values = c("Vol. 1 nd" = "#009E73",
               "Vol. 1 d" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines")))

# Print the plots
print(p_scat_comb.V1)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Scatterplot/prof_combined.Vol1Veff.png",
       plot = p_scat_comb.V1, width = 5, height = 5, dpi = 500)

# Vol 2
prof_combined.V2 <- prof.air.conc %>%
  select(congener, Conc.Air.V1.V2) %>%
  rename(Conc = Conc.Air.V1.V2) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V2.l) %>%
              rename(Conc = conc.V2.l) %>%
              mutate(Source = "Vol. 2 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V2.r) %>%
              rename(Conc = conc.V2.r) %>%
              mutate(Source = "Vol. 2 nd"))

# Plots
p_prof_comb.V2 <- ggplot(prof_combined.V2, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 2 nd" = "#009E73",
                               "Vol. 2 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(a)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V2)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal//Barplot/prof_combined.Vol2Veff.png",
       plot = p_prof_comb.V2, width = 10, height = 3, dpi = 500)

# Vol 3
prof_combined.V3 <- prof.air.conc %>%
  select(congener, Conc.Air.V3) %>%
  rename(Conc = Conc.Air.V3) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V3.l) %>%
              rename(Conc = conc.V3.l) %>%
              mutate(Source = "Vol. 3 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V3.r) %>%
              rename(Conc = conc.V3.r) %>%
              mutate(Source = "Vol. 3 d"))

# Plots
p_prof_comb.V3 <- ggplot(prof_combined.V3, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bar
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 3 nd" = "#009E73",
                               "Vol. 3 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(b)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V3)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol3Veff.png",
       plot = p_prof_comb.V3, width = 10, height = 3, dpi = 500)

# Vol 4
prof_combined.V4 <- prof.air.conc %>%
  select(congener, Conc.Air.V4) %>%
  rename(Conc = Conc.Air.V4) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V4.l) %>%
              rename(Conc = conc.V4.l) %>%
              mutate(Source = "Vol. 4 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V4.r) %>%
              rename(Conc = conc.V4.r) %>%
              mutate(Source = "Vol. 4 d"))

p_prof_comb.V4 <- ggplot(prof_combined.V4, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 4 nd" = "#009E73",
                               "Vol. 4 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(c)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V4)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol4Veff.png",
       plot = p_prof_comb.V4, width = 10, height = 3, dpi = 500)

# Vol 5
prof_combined.V5 <- prof.air.conc %>%
  select(congener, Conc.Air.V5) %>%
  rename(Conc = Conc.Air.V5) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V5.l) %>%
              rename(Conc = conc.V5.l) %>%
              mutate(Source = "Vol. 5 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V5.r) %>%
              rename(Conc = conc.V5.r) %>%
              mutate(Source = "Vol. 5 nd"))

p_prof_comb.V5 <- ggplot(prof_combined.V5, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 5 nd" = "#009E73",
                               "Vol. 5 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(d)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V5)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol5Veff.png",
       plot = p_prof_comb.V5, width = 10, height = 3, dpi = 500)

# Vol 6
prof_combined.V6 <- prof.air.conc %>%
  select(congener, Conc.Air.V6) %>%
  rename(Conc = Conc.Air.V6) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V6.l) %>%
              rename(Conc = conc.V6.l) %>%
              mutate(Source = "Vol. 6 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V6.r) %>%
              rename(Conc = conc.V6.r) %>%
              mutate(Source = "Vol. 6 d"))

# Plots
p_prof_comb.V6 <- ggplot(prof_combined.V6, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 6 nd" = "#009E73",
                               "Vol. 6 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(e)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V6)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol6Veff.png",
       plot = p_prof_comb.V6, width = 10, height = 3, dpi = 500)

# Vol 7
prof_combined.V7 <- prof.air.conc %>%
  select(congener, Conc.Air.V7) %>%
  rename(Conc = Conc.Air.V7) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V7.l) %>%
              rename(Conc = conc.V7.l) %>%
              mutate(Source = "Vol. 7 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V7.r) %>%
              rename(Conc = conc.V7.r) %>%
              mutate(Source = "Vol. 7 d"))

# Plots
p_prof_comb.V7 <- ggplot(prof_combined.V7, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0, 0.7), labels = scales::label_number(accuracy = 0.01)) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 7 nd" = "#009E73",
                               "Vol. 7 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(b)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V7)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol7VeffV2.png",
       plot = p_prof_comb.V7, width = 10, height = 3, dpi = 500)

# Convert to wide format
prof_wide <- prof_combined.V7 %>%
  pivot_wider(names_from = Source, values_from = Conc)

# Make a long-format for plotting multiple y-values against Air PCB
plot_data <- prof_wide %>%
  pivot_longer(
    cols = c(`Vol. 7 nd`, `Vol. 7 d`),  # use exact column names
    names_to = "Vol_Type",
    values_to = "Conc"
  )

ggplot(plot_data, aes(x = `Air PCB`, y = Conc, color = Vol_Type)) +
  geom_point(size = 2.5, shape = 21) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme_bw() +
  ylim(0, 0.5) +
  xlim(0, 0.5) +
  theme(
    aspect.ratio = 1,
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Air PCB Profile", y = "Volunteer PCB Profile",
       color = "Volunteer") +
  scale_color_manual(
    values = c("Vol. 7 nd" = "#009E73",
               "Vol. 7 d" = "#E69F00"),
    guide = guide_legend(key.size = unit(0.5, "lines"))
  )

# Vol 8
prof_combined.V8 <- prof.air.conc %>%
  select(congener, Conc.Air.V8) %>%
  rename(Conc = Conc.Air.V8) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, conc.V8.l) %>%
              rename(Conc = conc.V8.l) %>%
              mutate(Source = "Vol. 8 nd"))

# Plots
p_prof_comb.V8 <- ggplot(prof_combined.V8, aes(x = congener, y = Conc,
                                                  fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/20) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("Air PCB" = "blue",
                               "Vol. 8 nd" = "#009E73"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(1, 1),
        legend.justification = c(1 ,1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold")) +
  annotate("text", x = -Inf, y = Inf,
           label = "(f)", hjust = 0, vjust = 1, 
           size = 6, color = "black")

# Print the plots
print(p_prof_comb.V8)

# Save plot in folder
ggsave("Output/Plots/Profiles/Personal/Barplot/prof_combined.Vol8Veff.png",
       plot = p_prof_comb.V8, width = 10, height = 3, dpi = 500)

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

# Vol 4
prof_combined_wide.V4  <- prof_combined.V4 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V4[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V4 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V4

# Vol 5
prof_combined_wide.V5  <- prof_combined.V5 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V5[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V5 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V5

# Vol 6
prof_combined_wide.V6  <- prof_combined.V6 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V6[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V6 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V6

# Vol 7
prof_combined_wide.V7  <- prof_combined.V7 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V7[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V7 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V7

# Vol 8
prof_combined_wide.V8  <- prof_combined.V8 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V8[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V8 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V8

