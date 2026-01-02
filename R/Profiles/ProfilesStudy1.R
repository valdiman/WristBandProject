## Script to create PCB profiles &
# and calculate cosine theta for similarity analysis
# for Study 1

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
{
  PUF.0 <- read.csv("Data/IRO/SampleConcentrationStudy1.csv")
  PUF <- PUF.0[1:4, 8:180] # select rows and columns (ng/m3)
  WB.0 <- read.csv("Data/IRO/SampleMassStudy1.csv")
  WB <- WB.0[7:11, c(1, 5, 10:182)] # select rows and columns (ng/WB)
}

# Remove metadata
WB.1 <- subset(WB, select = -c(sid:time))

# Create PCB profiles -----------------------------------------------------
# (1) PUF
# Compute row sums
tmp.puf <- rowSums(PUF, na.rm = TRUE)
# Normalize each row so sum = 1
prof.puf.conc <- sweep(PUF, 1, tmp.puf, FUN = "/")
# Check
rowSums(prof.puf.conc, na.rm = TRUE)  # should be all 1

# (2) WBs
tmp.wb <- rowSums(WB.1, na.rm = TRUE)
prof.wb <- sweep(WB.1, 1, tmp.wb, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.wb, na.rm = TRUE)

# Plots -------------------------------------------------------------------
# PUF
df_long.puf <- prof.puf.conc %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(
    cols = matches("^PCB"),
    names_to = "PCB",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  mutate(
    PCB = factor(PCB, levels = unique(PCB)),
    Source = "Air PCB"
  )

df_summary <- df_long.puf %>%
  group_by(PCB) %>%
  summarise(
    mean_conc = mean(concentration),
    sd_conc   = sd(concentration),
    .groups = "drop"
  )

p.puf <- ggplot(df_summary, aes(x = PCB, y = mean_conc)) +
  geom_col(
    fill = "blue",
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(aes(ymin = mean_conc,
                    ymax = mean_conc + sd_conc),
                width = 0.2,
                color = "black") +
  theme_bw() +
  theme(
    aspect.ratio = 3/20,
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = NULL, y = expression(bold("Conc. Fraction " * Sigma * "PCB"))) +
  annotate("text", x = 150, y = 0.15, label = "Active low-volume sampler", 
    color = "black", size = 3, fontface = "bold")

p.puf

# Save plot in folder
ggsave("Output/Plots/Profiles/Study1/PUFConcentration.png", plot = p.puf,
       width = 10, height = 3, dpi = 500)

# WB
df_long.wb <- prof.wb %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(
    cols = matches("^PCB"),
    names_to = "PCB",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  mutate(
    PCB = factor(PCB, levels = unique(PCB)),
    Source = "WB PCB")

df_summary <- df_long.wb %>%
  group_by(PCB) %>%
  summarise(
    mean_conc = mean(concentration),
    sd_conc   = sd(concentration),
    .groups = "drop")

p.wb <- ggplot(df_summary, aes(x = PCB, y = mean_conc)) +
  geom_col( fill = "orange", color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean_conc, ymax = mean_conc + sd_conc),
                width = 0.2, color = "black") +
  theme_bw() +
  theme(
    aspect.ratio = 3/20,
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = NULL,y = expression(bold("Mass Fraction " * Sigma * "PCB"))) +
  annotate("text", x = 160, y = 0.15, label = "Static WB",  color = "black",
           size =3, fontface = "bold")

p.wb

# Save plot in folder
ggsave("Output/Plots/Profiles/Study1/WB.png", plot = p.wb,
       width = 10, height = 3, dpi = 500)

# Calculate cosine theta --------------------------------------------------
rownames(prof.puf.conc) <- sprintf("PUF%02d", 1:nrow(prof.puf.conc))
rownames(prof.wb) <- sprintf("WB%02d", 1:nrow(prof.wb))
combined <- bind_rows(prof.puf.conc, prof.wb)
combined_mat <- as.matrix(combined)
storage.mode(combined_mat) <- "numeric"
cosine_similarity <- cosine(t(combined_mat))
cosine_similarity[upper.tri(cosine_similarity, diag = TRUE)] <- NA
# View the resulting cosine similarity matrix
cosine_similarity
# Minimun cosine value
min(cosine_similarity, na.rm = TRUE)

