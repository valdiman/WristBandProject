## Script to fabricate PCB profiles &
# and calculate cosine theta for similarity analysis

# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("ggpubr")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(ggpubr) # ggarrange
}

# Sampling rates calculations under static conditions ---------------------
# Read data
{
  PUF <- read.csv("Data/PUF.csv") # ng/m3
  WB <- read.csv("Data/WB.csv")
}

# Remove metadata
PUF.1 <- subset(PUF, select = -c(sample))
WB.1 <- subset(WB, select = -c(sample:time))

# Create PCB profiles -----------------------------------------------------
# (1) PUF
# Compute row sums
row_totals <- rowSums(PUF.1, na.rm = TRUE)
# Normalize each row so sum = 1
prof.puf.conc <- sweep(PUF.1, 1, row_totals, FUN = "/")
# Check
rowSums(prof.puf.conc)  # should be all 1


# (2) WBs
tmp.wb <- rowSums(WB.1, na.rm = TRUE)
prof.wb <- sweep(WB.1, 1, tmp.wb, FUN = "/")
# Check sum of all PCBs (i.e., = 1)
rowSums(prof.wb, na.rm = TRUE)

# Plots -------------------------------------------------------------------
# PUF
df_long <- prof.puf.conc %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(
    cols = matches("^PCB"),   # matches any column starting with "PCB"
    names_to = "PCB",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  mutate(
    PCB = factor(PCB, levels = unique(PCB))  # keeps original order
  )


df_summary <- df_long %>%
  group_by(PCB) %>%
  summarise(
    mean_conc = mean(concentration),
    sd_conc   = sd(concentration),
    .groups = "drop"
  )

p.puf <- ggplot(df_summary, aes(x = PCB, y = mean_conc)) +
  geom_col(
    fill = "steelblue",
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(aes(ymin = mean_conc,  # start from mean
                    ymax = mean_conc + sd_conc),  # extend only upwards
                width = 0.2,
                color = "black") +
  theme_bw() +
  theme(
    aspect.ratio = 3/20,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Mean conc. fraction ",
  )

p.puf

# Save plot in folder
ggsave("Output/Plots/Profiles/Study1/PUFConcentration.png", plot = p.puf,
       width = 10, height = 3, dpi = 500)

# WB
df_long <- prof.wb %>%
  mutate(obs = 1:n()) %>%
  pivot_longer(
    cols = matches("^PCB"),   # matches any column starting with "PCB"
    names_to = "PCB",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  mutate(
    PCB = factor(PCB, levels = unique(PCB))  # keeps original order
  )


df_summary <- df_long %>%
  group_by(PCB) %>%
  summarise(
    mean_conc = mean(concentration),
    sd_conc   = sd(concentration),
    .groups = "drop"
  )

p.wb <- ggplot(df_summary, aes(x = PCB, y = mean_conc)) +
  geom_col(
    fill = "steelblue",
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(aes(ymin = mean_conc,  # start from mean
                    ymax = mean_conc + sd_conc),  # extend only upwards
                width = 0.2,
                color = "black") +
  theme_bw() +
  theme(
    aspect.ratio = 3/20,
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Mean conc. fraction ",
  )

p.wb

# Save plot in folder
ggsave("Output/Plots/Profiles/Study1/PUFConcentration.png", plot = p.puf,
       width = 10, height = 3, dpi = 500)
