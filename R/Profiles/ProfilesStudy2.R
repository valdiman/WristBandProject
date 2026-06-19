## Script to create PCB profiles &
# and calculate cosine theta for similarity analysis
# for Study 2 (mass)

# Install packages
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("lsa")
install.packages("SnowballC")

# Load libraries
{
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(lsa) # cosine theta function
}

# Read data ---------------------------------------------------------------
WB.s2 <- read.csv("Data/IRO/07_SampleWBMassStudy2.csv", check.names = FALSE)

# Volunteer 1 -------------------------------------------------------------
V1 <- WB.s2 %>%
  filter(volunteer == 1) %>%
  select(-c(sid, study, volunteer, wiped, vol.WB, area.WB))
V1.pcb <- subset(V1, select = -c(experiment:time))
# V1 static (n=3)
V1.static <- as.data.frame(t(colMeans(V1.pcb[1:3, ])))
# V1.nd last deployment 
V1.nd <- V1.pcb[8, ]
# V1.d last deployment
V1.d <- V1.pcb[13, ]
# Combine
V1.comb <- vctrs::vec_c(V1.static, V1.nd, V1.d)

# Create PCB profiles -----------------------------------------------------
tmp.V1 <- rowSums(V1.comb, na.rm = TRUE)
prof.V1 <- sweep(V1.comb, 1, tmp.V1, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.V1, na.rm = TRUE)
# Add row identifiers
prof.V1$Source <- c("Accumulated PCB Mass - Static WB",
                    "Accumulated PCB Mass - V1 nd",
                    "Accumulated PCB Mass - V1 d")

# Volunteer 1: Create PCB profile -----------------------------------------
pcb_order <- names(prof.V1)[names(prof.V1) != "Source"]

# Convert to long format
prof_combined.V1 <- prof.V1 %>%
  pivot_longer(
    cols = -Source,
    names_to = "congener",
    values_to = "Conc") %>%
  mutate(
    congener = factor(congener, levels = pcb_order))

prof_combined.V1$congener <- factor(
  prof_combined.V1$congener,
  levels = pcb_order)

# Volunteer 1: Create plot ------------------------------------------------
p_prof_comb.V1 <- ggplot(prof_combined.V1, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black",
           linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Mass Fraction " * Sigma * "PCB"))) +
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
  scale_fill_manual(values = c("Accumulated PCB Mass - Static WB" = "blue",
                               "Accumulated PCB Mass - V1 nd" = "#009E73",
                               "Accumulated PCB Mass - V1 d" = "#E69F00")) +
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See plot
p_prof_comb.V1

# Save plot in folder
ggsave("Output/Plots/Profiles/Study2/V1_study2.png", plot = p_prof_comb.V1,
       width = 22, height = 5, dpi = 500)

# Volunteer 1: Calculate cosine theta -------------------------------------
prof_combined_wide.V1  <- prof_combined.V1 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V1[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V1 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V1

# Volunteer 2 -------------------------------------------------------------
V2 <- WB.s2 %>%
  filter(volunteer == 2) %>%
  select(-c(sid, study, volunteer, wiped, vol.WB, area.WB))
V2.pcb <- subset(V2, select = -c(experiment:time))
# V1 static (n=3)
V2.static <- as.data.frame(t(colMeans(V2.pcb[1:3, ])))
# V2.d last deployment
V2.d <- V2.pcb[8, ]
# Combine
V2.comb <- vctrs::vec_c(V2.static, V2.d)

# Create PCB profiles -----------------------------------------------------
tmp.V2 <- rowSums(V2.comb, na.rm = TRUE)
prof.V2 <- sweep(V2.comb, 1, tmp.V2, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.V2, na.rm = TRUE)
# Add row identifiers
prof.V2$Source <- c("Accumulated PCB Mass - Static WB", "Accumulated PCB Mass - V2 d")

# Volunteer 2: Create PCB profile -----------------------------------------
pcb_order <- names(prof.V2)[names(prof.V2) != "Source"]

# Convert to long format
prof_combined.V2 <- prof.V2 %>%
  pivot_longer(
    cols = -Source,
    names_to = "congener",
    values_to = "Conc") %>%
  mutate(
    congener = factor(congener, levels = pcb_order))

prof_combined.V2$congener <- factor(
  prof_combined.V2$congener,
  levels = pcb_order)

# Volunteer 2: Create plot ------------------------------------------------
p_prof_comb.V2 <- ggplot(prof_combined.V2, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black",
           linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Mass Fraction " * Sigma * "PCB"))) +
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
  scale_fill_manual(values = c("Accumulated PCB Mass - Static WB" = "blue",
                               "Accumulated PCB Mass - V2 d" = "#E69F00")) +
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See plot
p_prof_comb.V2

# Save plot in folder
ggsave("Output/Plots/Profiles/Study2/V2_study2.png", plot = p_prof_comb.V2,
       width = 22, height = 5, dpi = 500)

# Volunteer 2: Calculate cosine theta -------------------------------------
prof_combined_wide.V2  <- prof_combined.V2 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V2[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V2 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V2

# Volunteer 3 1st week-------------------------------------------------------------
V3 <- WB.s2 %>%
  filter(volunteer == 3) %>%
  select(-c(sid, study, volunteer, wiped, vol.WB, area.WB))
V3.pcb <- subset(V3, select = -c(experiment:time))
# V3 static last deployment 1st week
V3.static <- V3.pcb[3, ]
# V3.nd last deployment 1st week 
V3.nd.1w <- V3.pcb[9, ]
# Combine
V3.comb <- vctrs::vec_c(V3.static, V3.nd.1w)

# Create PCB profiles -----------------------------------------------------
tmp.V3 <- rowSums(V3.comb, na.rm = TRUE)
prof.V3 <- sweep(V3.comb, 1, tmp.V3, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.V3, na.rm = TRUE)
# Add row identifiers
prof.V3$Source <- c("Accumulated PCB Mass - Static WB 1st week",
                    "Accumulated PCB Mass - V3 nd 1st week")

# Volunteer 3 1st week: Create PCB profile -----------------------------------------
pcb_order <- names(prof.V3)[names(prof.V3) != "Source"]

# Convert to long format
prof_combined.V3 <- prof.V3 %>%
  pivot_longer(
    cols = -Source,
    names_to = "congener",
    values_to = "Conc") %>%
  mutate(
    congener = factor(congener, levels = pcb_order))

prof_combined.V3$congener <- factor(
  prof_combined.V3$congener,
  levels = pcb_order)

# Volunteer 3 1st week: Create plot ------------------------------------------------
p_prof_comb.V3 <- ggplot(prof_combined.V3, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black",
           linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Mass Fraction " * Sigma * "PCB"))) +
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
  scale_fill_manual(values = c("Accumulated PCB Mass - Static WB 1st week" = "blue",
                               "Accumulated PCB Mass - V3 nd 1st week" = "#009E73")) +
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See plot
p_prof_comb.V3

# Save plot in folder
ggsave("Output/Plots/Profiles/Study2/V3_1stW_study2.png", plot = p_prof_comb.V3,
       width = 22, height = 5, dpi = 500)

# Volunteer 3 1st week: Calculate cosine theta -------------------------------------
prof_combined_wide.V3  <- prof_combined.V3 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V3[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V3 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V3

# Volunteer 3 2nd week-------------------------------------------------------------
# V3 static last deployment 2nd week
V3.static <- V3.pcb[6, ]
# V3.nd last deployment 2nd week 
V3.nd.2w <- V3.pcb[12, ]
# Combine
V3.comb <- vctrs::vec_c(V3.static, V3.nd.2w)

# Create PCB profiles -----------------------------------------------------
tmp.V3 <- rowSums(V3.comb, na.rm = TRUE)
prof.V3 <- sweep(V3.comb, 1, tmp.V3, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.V3, na.rm = TRUE)
# Add row identifiers
prof.V3$Source <- c("Accumulated PCB Mass - Static WB 2nd week",
                    "Accumulated PCB Mass - V3 nd 2nd week")

# Volunteer 3 2nd week: Create PCB profile -----------------------------------------
pcb_order <- names(prof.V3)[names(prof.V3) != "Source"]

# Convert to long format
prof_combined.V3 <- prof.V3 %>%
  pivot_longer(
    cols = -Source,
    names_to = "congener",
    values_to = "Conc") %>%
  mutate(
    congener = factor(congener, levels = pcb_order))

prof_combined.V3$congener <- factor(
  prof_combined.V3$congener,
  levels = pcb_order)

# Volunteer 3 2nd week: Create plot ------------------------------------------------
p_prof_comb.V3 <- ggplot(prof_combined.V3, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black",
           linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Mass Fraction " * Sigma * "PCB"))) +
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
  scale_fill_manual(values = c("Accumulated PCB Mass - Static WB 2nd week" = "blue",
                               "Accumulated PCB Mass - V3 nd 2nd week" = "#009E73")) +
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See plot
p_prof_comb.V3

# Save plot in folder
ggsave("Output/Plots/Profiles/Study2/V3_2ndW_study2.png", plot = p_prof_comb.V3,
       width = 22, height = 5, dpi = 500)

# Volunteer 3 2nd week: Calculate cosine theta -------------------------------------
prof_combined_wide.V3  <- prof_combined.V3 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V3[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V3 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V3

# Volunteer 3 non-wiped & wiped-------------------------------------------------------------
# V3 static last deployment no-wiped
V3.static <- V3.pcb[15, ]
# V3.d last deployment no-wiped
V3.d.2nw <- V3.pcb[18, ]
# V3.d last deployment wiped
V3.d.2w <- V3.pcb[21, ]
# Combine
V3.comb <- vctrs::vec_c(V3.static, V3.d.2nw, V3.d.2w)

# Create PCB profiles -----------------------------------------------------
tmp.V3 <- rowSums(V3.comb, na.rm = TRUE)
prof.V3 <- sweep(V3.comb, 1, tmp.V3, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.V3, na.rm = TRUE)
# Add row identifiers
prof.V3$Source <- c("Accumulated PCB Mass - Static WB",
                    "Accumulated PCB Mass - V3 d Non-wiped",
                    "Accumulated PCB Mass - V3 d Wiped")

# Volunteer 3 non-wiped & wiped: Create PCB profile -----------------------------------------
pcb_order <- names(prof.V3)[names(prof.V3) != "Source"]

# Convert to long format
prof_combined.V3 <- prof.V3 %>%
  pivot_longer(
    cols = -Source,
    names_to = "congener",
    values_to = "Conc") %>%
  mutate(
    congener = factor(congener, levels = pcb_order))

prof_combined.V3$congener <- factor(
  prof_combined.V3$congener,
  levels = pcb_order)

# Volunteer 3 non-wiped & wiped: Create plot ------------------------------------------------
p_prof_comb.V3 <- ggplot(prof_combined.V3, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(stat = "identity", width = 1, color = "black",
           linewidth = 0.2) +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(limits = c(0, 0.15), n.breaks = 3) +
  theme_bw() +
  ylab(expression(bold("Mass Fraction " * Sigma * "PCB"))) +
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
  scale_fill_manual(values = c("Accumulated PCB Mass - Static WB" = "blue",
                               "Accumulated PCB Mass - V3 d Non-wiped" = "#009E73",
                               "Accumulated PCB Mass - V3 d Wiped" = "#E69F00")) +
  scale_x_discrete(
    labels = function(x) gsub("\\.", "+", x))

# See plot
p_prof_comb.V3

# Save plot in folder
ggsave("Output/Plots/Profiles/Study2/V3_wipe_study2.png", plot = p_prof_comb.V3,
       width = 22, height = 5, dpi = 500)

# Volunteer 3 non-wiped & wiped: Calculate cosine theta -------------------------------------
prof_combined_wide.V3  <- prof_combined.V3 %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.V3[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.V3 <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.V3
