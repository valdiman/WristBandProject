# Concentration estimation

# Install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("scales")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("ggpubr")
install.packages("lsa")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(ggpubr) # ggarrange
  library(lsa) # cosine theta function
}

# Read measured values from excel -----------------------------------------
data <- data.frame(read_excel("Data/Volunteers.xlsx", sheet = "Sheet1",
                                       col_names = TRUE, col_types = NULL))

# Mass accumulated in static WBs ------------------------------------------
{
  # For MI & YA
  mass.Mi.Ya.Stat <- data[8, 4:176]
  # For EA
  mass.Ea.Stat <- data[7, 4:176]
  # For Cr
  mass.Cr.Stat <- data.frame(t(colMeans(data[14:15, 4:176])))
  # For Hu
  mass.Hu.Stat <- data[13, 4:176]
  # For Xu
  mass.Xu.Stat <- data.frame(t(colMeans(data[19:21, 4:176])))
  # Combine concentrations
  wb.stat <- rbind(mass.Mi.Ya.Stat, mass.Ea.Stat, mass.Cr.Stat, mass.Hu.Stat, mass.Xu.Stat)
  # Add row names
  rownames(wb.stat) <- c("Mass.Air.Mi.Ya.Stat", "Mass.Air.Ea.Stat", "Mass.Air.Cr.Stat",
                          "Mass.Air.Hu.Stat", "Mass.Air.Xu.Stat")
}

# Mass accumulated in wore WBs --------------------------------------------
{
  wb.Mi.l <- data[1, 4:176]
  wb.Mi.r <- data[2, 4:176]
  wb.Ya.l <- data[3, 4:176]
  wb.Ya.r <- data[4, 4:176]
  wb.Ea.l <- data[5, 4:176]
  wb.Ea.r <- data[6, 4:176]
  wb.Cr.l <- data[9, 4:176]
  wb.Cr.r <- data[10, 4:176]
  wb.Hu.l <- data[11, 4:176]
  wb.Hu.r <- data[12, 4:176]
  wb.Xu.l <- data[17, 4:176]
  wb.Xu.r <- data[18, 4:176]
  # Combined wore WBs
  wb.wr <- rbind(wb.Mi.l, wb.Mi.r, wb.Ya.l, wb.Ya.r, wb.Ea.l, wb.Ea.r,
                 wb.Cr.l, wb.Cr.r, wb.Hu.l, wb.Hu.r, wb.Xu.l, wb.Xu.r)
  # Add row names
  rownames(wb.wr) <- c("Mass.Mi.l.wr", "Mass.Mi.r.wr", "Mass.Ya.l.wr", "Mass.Ya.r.wr",
                       "Mass.Ea.l.wr", "Mass.Ea.r.wr", "Mass.Cr.l.wr", "Mass.Cr.r.wr",
                       "Mass.Hu.l.wr", "Mass.Hu.r.wr", "Mass.Xu.l.wr", "Mass.Xu.r.wr")
}

# Create mass PCB congener profiles ---------------------------------------
# (1) Air WBs
tmp.wb.stat <- rowSums(wb.stat, na.rm = TRUE)
prof.wb.stat <- sweep(wb.stat, 1, tmp.wb.stat, FUN = "/")
prof.wb.stat <- data.frame(t(prof.wb.stat))
congener <- rownames(prof.wb.stat)
prof.wb.stat <- cbind(congener, prof.wb.stat)
rownames(prof.wb.stat) <- NULL
prof.wb.stat[, 2:6] <- lapply(prof.wb.stat[, 2:6], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.stat$congener <- factor(prof.wb.stat$congener,
                                levels = unique(prof.wb.stat$congener))

# (2) Wore WBs
tmp.wb.wr <- rowSums(wb.wr, na.rm = TRUE)
prof.wb.wr <- sweep(wb.wr, 1, tmp.wb.wr, FUN = "/")
prof.wb.wr <- data.frame(t(prof.wb.wr))
congener <- rownames(prof.wb.wr)
prof.wb.wr <- cbind(congener, prof.wb.wr)
rownames(prof.wb.wr) <- NULL
prof.wb.wr[, 2:13] <- lapply(prof.wb.wr[, 2:13], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.wr$congener <- factor(prof.wb.wr$congener,
                                levels = unique(prof.wb.wr$congener))

# Mass profile plots ------------------------------------------------------
# Mi
prof_combined.Mi <- prof.wb.stat %>%
  select(congener, Mass.Air.Mi.Ya.Stat) %>%
  rename(Mass = Mass.Air.Mi.Ya.Stat) %>%
  mutate(Source = "prof.wb.stat.Mi") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Mi.l.wr) %>%
              rename(Mass = Mass.Mi.l.wr) %>%
              mutate(Source = "prof.wb.wr.Mi.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Mi.r.wr) %>%
              rename(Mass = Mass.Mi.r.wr) %>%
              mutate(Source = "prof.wb.wr.Mi.r"))

# Plots
ggplot(prof_combined.Mi, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Mi" = "#0072B2",  # Blue
                               "prof.wb.wr.Mi.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Mi.r" = "#009E73")) # Green

# Ya
prof_combined.Ya <- prof.wb.stat %>%
  select(congener, Mass.Air.Mi.Ya.Stat) %>%
  rename(Mass = Mass.Air.Mi.Ya.Stat) %>%
  mutate(Source = "prof.wb.stat.Ya") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ya.l.wr) %>%
              rename(Mass = Mass.Ya.l.wr) %>%
              mutate(Source = "prof.wb.wr.Ya.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ya.r.wr) %>%
              rename(Mass = Mass.Ya.r.wr) %>%
              mutate(Source = "prof.wb.wr.Ya.r"))

# Plots
ggplot(prof_combined.Ya, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Ya" = "#0072B2",  # Blue
                               "prof.wb.wr.Ya.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Ya.r" = "#009E73")) # Green

# Ea
prof_combined.Ea <- prof.wb.stat %>%
  select(congener, Mass.Air.Ea.Stat) %>%
  rename(Mass = Mass.Air.Ea.Stat) %>%
  mutate(Source = "prof.wb.stat.Ea") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ea.l.wr) %>%
              rename(Mass = Mass.Ea.l.wr) %>%
              mutate(Source = "prof.wb.wr.Ea.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ea.r.wr) %>%
              rename(Mass = Mass.Ea.r.wr) %>%
              mutate(Source = "prof.wb.wr.Ea.r"))

# Plots
ggplot(prof_combined.Ea, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Ea" = "#0072B2",  # Blue
                               "prof.wb.wr.Ea.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Ea.r" = "#009E73")) # Green

# Cr
prof_combined.Cr <- prof.wb.stat %>%
  select(congener, Mass.Air.Cr.Stat) %>%
  rename(Mass = Mass.Air.Cr.Stat) %>%
  mutate(Source = "prof.wb.stat.Cr") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Cr.l.wr) %>%
              rename(Mass = Mass.Cr.l.wr) %>%
              mutate(Source = "prof.wb.wr.Cr.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Cr.r.wr) %>%
              rename(Mass = Mass.Cr.r.wr) %>%
              mutate(Source = "prof.wb.wr.Cr.r"))

# Plots
ggplot(prof_combined.Cr, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Cr" = "#0072B2",  # Blue
                               "prof.wb.wr.Cr.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Cr.r" = "#009E73")) # Green

# Hu
prof_combined.Hu <- prof.wb.stat %>%
  select(congener, Mass.Air.Hu.Stat) %>%
  rename(Mass = Mass.Air.Hu.Stat) %>%
  mutate(Source = "prof.wb.stat.Hu") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Hu.l.wr) %>%
              rename(Mass = Mass.Hu.l.wr) %>%
              mutate(Source = "prof.wb.wr.Hu.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Hu.r.wr) %>%
              rename(Mass = Mass.Hu.r.wr) %>%
              mutate(Source = "prof.wb.wr.Hu.r"))

# Plots
ggplot(prof_combined.Hu, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Hu" = "#0072B2",  # Blue
                               "prof.wb.wr.Hu.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Hu.r" = "#009E73")) # Green

# Xu
prof_combined.Xu <- prof.wb.stat %>%
  select(congener, Mass.Air.Xu.Stat) %>%
  rename(Mass = Mass.Air.Xu.Stat) %>%
  mutate(Source = "prof.wb.stat.Xu") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Xu.l.wr) %>%
              rename(Mass = Mass.Xu.l.wr) %>%
              mutate(Source = "prof.wb.wr.Xu.l")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Xu.r.wr) %>%
              rename(Mass = Mass.Xu.r.wr) %>%
              mutate(Source = "prof.wb.wr.Xu.r"))

# Plots
pXu <- ggplot(prof_combined.Xu, aes(x = congener, y = Mass, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1) +
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("prof.wb.stat.Xu" = "#0072B2",  # Blue
                               "prof.wb.wr.Xu.l" = "#E69F00",  # Orange
                               "prof.wb.wr.Xu.r" = "#009E73")) # Green

# Print the plot
print(pXu)

# Save plot in folder
ggsave("Output/Plots/Profiles/MassProfXu.png", plot = pXu, width = 6,
       height = 4, dpi = 500)

# Create concentration PCB congener profiles ------------------------------
# (1) Air WBs
# = massWB/(0.5*time.day)
time <- data$time.day[c(8, 7, 14, 13, 19)]
conc.wb.air <- as.data.frame(wb.stat/(0.5*time))
rownames(conc.wb.air) <- c("Conc.Air.Mi.Ya.Stat",
                                "Conc.Air.Ea.Stat", "Conc.Air.Cr.Stat",
                                "Conc.Air.Hu.Stat", "Conc.Air.Xu.Stat")
# Check total PCB
tPCB.conc.air <- rowSums(conc.wb.air, na.rm = TRUE)
# See
tPCB.conc.air

# (2) Wore WBs
# Read calculated average sampling rates
sr <- read.csv("Output/Data/csv/Ave.SRs.csv")
# Select only average sampling rate
sr <- sr[, 1:2]

# Extract congener names from sr
congener_names <- sr$congener
# Subset wb.wr to include only the common congeners
wb_common <- wb.wr[, intersect(colnames(wb.wr), congener_names)]
# Subset sr to include only the common congeners
sr_common <- sr[sr$congener %in% colnames(wb_common), ]
# Ensure that sr_common is ordered to match the columns of wb_common
sr_common <- sr_common[match(colnames(wb_common), sr_common$congener), ]
# Divide each element in wb_common by corresponding element in sr_common$Average_Sampling_Rate
wb_div_sr <- sweep(wb_common, 2, sr_common$Average_Sampling_Rate, FUN = "/")
# !!! Use a constant sampling rate
# wb_div_sr <- sweep(wb_common, 2, 1.5, FUN = "/")
# Extract time.day from wb.wr
wb_time_day <- data$time.day[c(1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 17, 18)]
# Divide wb_div_sr further by corresponding value in wb_time_day
conc.wb.wr <- sweep(wb_div_sr, 1, wb_time_day, FUN = "/")
rownames(conc.wb.wr) <- c("Conc.Mi.l.wr", "Conc.Mi.r.wr", "Conc.Ya.l.wr",
                          "Conc.Ya.r.wr", "Conc.Ea.l.wr",
                               "Conc.Ea.r.wr", "Conc.Cr.l.wr", "Conc.Cr.r.wr",
                               "Conc.Hu.l.wr", "Conc.Hu.r.wr", "Conc.Xu.l.wr",
                               "Conc.Xu.r.wr")
# Check total PCB
tPCB.conc.wb.wr <- rowSums(conc.wb.wr, na.rm = TRUE)
# See
tPCB.conc.wb.wr

# Transpose conc.wb.air
conc.wb.air.t <- data.frame(t(conc.wb.air))
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb.wr), rownames(conc.wb.air.t))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.wb.air.t))
# Subset conc.air to include only the rows with matching row names
conc.wb.air.t <- conc.wb.air.t[matching_indices, ]
# Transpose conc.wb.air.t again
conc.wb.air.t <- data.frame(t(conc.wb.air.t))

# Profiles
# (1) Air WBs
tmp.wb.air <- rowSums(conc.wb.air.t, na.rm = TRUE)
prof.wb.air.conc <- sweep(conc.wb.air.t, 1, tmp.wb.air, FUN = "/")
prof.wb.air.conc <- data.frame(t(prof.wb.air.conc))
congener <- rownames(prof.wb.air.conc)
prof.wb.air.conc <- cbind(congener, prof.wb.air.conc)
rownames(prof.wb.air.conc) <- NULL
prof.wb.air.conc[, 2:6] <- lapply(prof.wb.air.conc[, 2:6], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.air.conc$congener <- factor(prof.wb.air.conc$congener,
                                levels = unique(prof.wb.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.air.conc[, 2:6], na.rm = TRUE)

# (2) Wore WBs
tmp.wb.wr <- rowSums(conc.wb.wr, na.rm = TRUE)
prof.wb.wr.conc <- sweep(conc.wb.wr, 1, tmp.wb.wr, FUN = "/")
prof.wb.wr.conc <- data.frame(t(prof.wb.wr.conc))
congener <- rownames(prof.wb.wr.conc)
prof.wb.wr.conc <- cbind(congener, prof.wb.wr.conc)
rownames(prof.wb.wr.conc) <- NULL
prof.wb.wr.conc[, 2:13] <- lapply(prof.wb.wr.conc[, 2:13], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.wr.conc$congener <- factor(prof.wb.wr.conc$congener,
                               levels = unique(prof.wb.wr.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.wr.conc[, 2:13], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# Mi (1)
prof_combined.Mi <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Mi.Ya.Stat) %>%
  rename(Conc = Conc.Air.Mi.Ya.Stat) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Mi.l.wr) %>%
              rename(Conc = Conc.Mi.l.wr) %>%
              mutate(Source = "Vol. 1 nd")) %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Mi.r.wr) %>%
              rename(Conc = Conc.Mi.r.wr) %>%
              mutate(Source = "Vol. 1 d"))  # Change this label

# Plot
# Create the plot with the legend moved inside
prof_combined.Mi <- ggplot(prof_combined.Mi, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           size = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Mi)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol1.png", plot = prof_combined.Mi,
       width = 10, height = 5, dpi = 500)

# Ya
prof_combined.Ya <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Mi.Ya.Stat) %>%
  rename(Conc = Conc.Air.Mi.Ya.Stat) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Ya.l.wr) %>%
              rename(Conc = Conc.Ya.l.wr) %>%
              mutate(Source = "Vol. 2 d")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Ya.r.wr) %>%
              rename(Conc = Conc.Ya.r.wr) %>%
              mutate(Source = "Vol. 2 nd"))

# Plots
prof_combined.Ya <- ggplot(prof_combined.Ya, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           size = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Ya)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol2.png", plot = prof_combined.Ya,
       width = 10, height = 5, dpi = 500)

# Ea
prof_combined.Ea <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Ea.Stat) %>%
  rename(Conc = Conc.Air.Ea.Stat) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Ea.l.wr) %>%
              rename(Conc = Conc.Ea.l.wr) %>%
              mutate(Source = "Vol. 3 nd")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Ea.r.wr) %>%
              rename(Conc = Conc.Ea.r.wr) %>%
              mutate(Source = "Vol. 3 d"))

# Plots
prof_combined.Ea <- ggplot(prof_combined.Ea, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the barshttp://127.0.0.1:8373/graphics/plot_zoom_png?width=1872&height=900
           size = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Ea)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol3.png", plot = prof_combined.Ea,
       width = 10, height = 5, dpi = 500)

# Cr
prof_combined.Cr <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Cr.Stat) %>%
  rename(Conc = Conc.Air.Cr.Stat) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Cr.l.wr) %>%
              rename(Conc = Conc.Cr.l.wr) %>%
              mutate(Source = "Vol. 4 nd")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Cr.r.wr) %>%
              rename(Conc = Conc.Cr.r.wr) %>%
              mutate(Source = "Vol. 4 d"))

prof_combined.Cr <- ggplot(prof_combined.Cr, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           size = 0.2) +  # Set the thickness of the black edges (fine line)http://127.0.0.1:8373/graphics/plot_zoom_png?width=1856&height=861
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Cr)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol4.png", plot = prof_combined.Cr,
       width = 10, height = 5, dpi = 500)

# Hu
prof_combined.Hu <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Hu.Stat) %>%
  rename(Conc = Conc.Air.Hu.Stat) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Hu.l.wr) %>%
              rename(Conc = Conc.Hu.l.wr) %>%
              mutate(Source = "Vol. 5 d")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Hu.r.wr) %>%
              rename(Conc = Conc.Hu.r.wr) %>%
              mutate(Source = "Vol. 5 nd"))

prof_combined.Hu <- ggplot(prof_combined.Hu, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           size = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Hu)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol5.png", plot = prof_combined.Hu,
       width = 10, height = 5, dpi = 500)

# Xu
prof_combined.Xu <- prof.wb.air.conc %>%
  select(congener, Conc.Air.Xu.Stat) %>%
  rename(Conc = Conc.Air.Xu.Stat) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Xu.l.wr) %>%
              rename(Conc = Conc.Xu.l.wr) %>%
              mutate(Source = "Vol. 6 nd")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, Conc.Xu.r.wr) %>%
              rename(Conc = Conc.Xu.r.wr) %>%
              mutate(Source = "Vol. 6 d"))

# Plots
prof_combined.Xu <- ggplot(prof_combined.Xu, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           size = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.15) +
  theme_bw() +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Concentration fraction "*Sigma*"PCB"))) +
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
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(prof_combined.Xu)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol6.png", plot = prof_combined.Xu,
       width = 10, height = 5, dpi = 500)

# Calculate cosine theta --------------------------------------------------
# Need to change the format of prof_combined...
# Create a term-document matrix from the 'Source' variable
# Mi
prof_combined_wide.Mi  <- prof_combined.Mi %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Mi[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Mi <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Mi

# Ya
prof_combined_wide.Ya  <- prof_combined.Ya %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Ya[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Ya <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Ya

# Ea
prof_combined_wide.Ea  <- prof_combined.Ea %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Ea[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Ea <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Ea

# Cr
prof_combined_wide.Cr  <- prof_combined.Cr %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Cr[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Cr <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Cr

# Hu
prof_combined_wide.Hu  <- prof_combined.Hu %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Hu[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Hu <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Hu

# Xu
prof_combined_wide.Xu  <- prof_combined.Xu %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Xu[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Xu <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Xu

