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
data <- data.frame(read_excel("Data/VolunteersV02.xlsx", sheet = "Sheet1",
                              col_names = TRUE, col_types = NULL))

# Calculate air PCB concentration from static WBs -------------------------
{
  # Stat samples
  # 1 and 2 are from the same space
  # 3 is for An and Xu
  mass.stat.12 <- data.frame(t(colMeans(data[1:2, 3:175])))
  mass.stat.3 <- data[3, 3:175]
  # Bind the data frames together row-wise
  wb.stat <- rbind(mass.stat.12, mass.stat.3)
  # Change column names of the  columns
  rownames(wb.stat) <- c("Mass.stat.1", "Mass.stat.2")
}

# Mass accumulated in wore WBs --------------------------------------------
{
  wb.Mi.o <- data[4, 3:175]
  wb.Mi.h <- data[5, 3:175]
  wb.Ea.h <- data[6, 3:175]
  wb.Ea.o <- data[7, 3:175]
  wb.Ya.o <- data[8, 3:175]
  wb.Ya.h <- data[9, 3:175]
  wb.An.h <- data[10, 3:175]
  wb.An.o <- data[11, 3:175]
  wb.Xu.h <- data[12, 3:175]
  wb.Xu.o <- data[13, 3:175]
 # Combined wore WBs
 wb.wr <- rbind(wb.Mi.o, wb.Mi.h, wb.Ea.o, wb.Ea.h, wb.Ya.o, wb.Ya.h,
               wb.An.o, wb.An.h, wb.Xu.o, wb.Xu.h)
 # Add row names
 rownames(wb.wr) <- c("Mass.Mi.o", "Mass.Mi.h", "Mass.Ea.o", "Mass.Ea.h",
                       "Mass.Ya.o", "Mass.Ya.h", "Mass.An.o", "Mass.An.h",
                       "Mass.Xu.o", "Mass.Xu.h")
}

# Create mass PCB congener profiles ---------------------------------------
# (1) Air WBs
tmp.wb.stat <- rowSums(wb.stat, na.rm = TRUE)
prof.wb.stat <- sweep(wb.stat, 1, tmp.wb.stat, FUN = "/")
prof.wb.stat <- data.frame(t(prof.wb.stat))
congener <- rownames(prof.wb.stat)
prof.wb.stat <- cbind(congener, prof.wb.stat)
rownames(prof.wb.stat) <- NULL
prof.wb.stat[, 2:3] <- lapply(prof.wb.stat[, 2:3], as.numeric)
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
prof.wb.wr[, 2:10] <- lapply(prof.wb.wr[, 2:10], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.wr$congener <- factor(prof.wb.wr$congener,
                                levels = unique(prof.wb.wr$congener))

# Mass profile plots ------------------------------------------------------
# Mi
prof_combined.Mi <- prof.wb.stat %>%
  select(congener, Mass.stat.1) %>%
  rename(Mass = Mass.stat.1) %>%
  mutate(Source = "prof.wb.stat.Mi") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Mi.o) %>%
              rename(Mass = Mass.Mi.o) %>%
              mutate(Source = "prof.wb.wr.Mi.o")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Mi.h) %>%
              rename(Mass = Mass.Mi.h) %>%
              mutate(Source = "prof.wb.wr.Mi.h"))

# Plots
ggplot(prof_combined.Mi, aes(x = congener, y = Mass, fill = Source)) +
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
  scale_fill_manual(values = c("prof.wb.stat.Mi" = "#0072B2",  # Blue
                               "prof.wb.wr.Mi.o" = "#E69F00",  # Orange
                               "prof.wb.wr.Mi.h" = "#009E73")) # Green

# Ea
prof_combined.Ea <- prof.wb.stat %>%
  select(congener, Mass.stat.1) %>%
  rename(Mass = Mass.stat.1) %>%
  mutate(Source = "prof.wb.stat.Ea") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ea.o) %>%
              rename(Mass = Mass.Ea.o) %>%
              mutate(Source = "prof.wb.wr.Ea.o")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ea.h) %>%
              rename(Mass = Mass.Ea.h) %>%
              mutate(Source = "prof.wb.wr.Ea.h"))

# Plots
ggplot(prof_combined.Ea, aes(x = congener, y = Mass, fill = Source)) +
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
  scale_fill_manual(values = c("prof.wb.stat.Ea" = "#0072B2",  # Blue
                               "prof.wb.wr.Ea.o" = "#E69F00",  # Orange
                               "prof.wb.wr.Ea.h" = "#009E73")) # Green

# Ya
prof_combined.Ya <- prof.wb.stat %>%
  select(congener, Mass.stat.1) %>%
  rename(Mass = Mass.stat.1) %>%
  mutate(Source = "prof.wb.stat.Ya") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ya.o) %>%
              rename(Mass = Mass.Ya.o) %>%
              mutate(Source = "prof.wb.wr.Ya.o")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Ya.h) %>%
              rename(Mass = Mass.Ya.h) %>%
              mutate(Source = "prof.wb.wr.Ya.h"))

# Plots
ggplot(prof_combined.Ya, aes(x = congener, y = Mass, fill = Source)) +
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
  scale_fill_manual(values = c("prof.wb.stat.Ya" = "#0072B2",  # Blue
                               "prof.wb.wr.Ya.o" = "#E69F00",  # Orange
                               "prof.wb.wr.Ya.h" = "#009E73")) # Green

# An
prof_combined.An <- prof.wb.stat %>%
  select(congener, Mass.stat.2) %>%
  rename(Mass = Mass.stat.2) %>%
  mutate(Source = "prof.wb.stat.An") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.An.o) %>%
              rename(Mass = Mass.An.o) %>%
              mutate(Source = "prof.wb.wr.An.o")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.An.h) %>%
              rename(Mass = Mass.An.h) %>%
              mutate(Source = "prof.wb.wr.An.h"))

# Plots
ggplot(prof_combined.An, aes(x = congener, y = Mass, fill = Source)) +
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
  scale_fill_manual(values = c("prof.wb.stat.An" = "#0072B2",  # Blue
                               "prof.wb.wr.An.o" = "#E69F00",  # Orange
                               "prof.wb.wr.An.h" = "#009E73")) # Green

# Xu
prof_combined.Xu <- prof.wb.stat %>%
  select(congener, Mass.stat.2) %>%
  rename(Mass = Mass.stat.2) %>%
  mutate(Source = "prof.wb.stat.Xu") %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Xu.o) %>%
              rename(Mass = Mass.Xu.o) %>%
              mutate(Source = "prof.wb.wr.Xu.o")) %>%
  bind_rows(prof.wb.wr %>%
              select(congener, Mass.Xu.h) %>%
              rename(Mass = Mass.Xu.h) %>%
              mutate(Source = "prof.wb.wr.Xu.h"))

# Plots
pXu <- ggplot(prof_combined.Xu, aes(x = congener, y = Mass, fill = Source)) +
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
  scale_fill_manual(values = c("prof.wb.stat.Xu" = "#0072B2",  # Blue
                               "prof.wb.wr.Xu.o" = "#E69F00",  # Orange
                               "prof.wb.wr.Xu.h" = "#009E73")) # Green

# Print the plot
print(pXu)

# Save plot in folder for example
ggsave("Output/Plots/Profiles/OfficeHome/MassProfXu.png", plot = pXu, width = 6,
       height = 4, dpi = 500)

# Estimate air concentrations ---------------------------------------------
# (1) Air WBs
# = massWB/(0.5*time.day)
data.Stat.1 <- data[1, 3:175]
data.Stat.2 <- data[2, 3:175]
data.Stat.3 <- data[3, 3:175]
# Calculate air concentration in ng/m3
# = massWB/(0.5*time.day)
# i and ii are from the same office
conc.air.i <- as.data.frame(data.Stat.1/(0.5*data$time.day[1]))
conc.air.ii <- as.data.frame(data.Stat.2/(0.5*data$time.day[2]))
conc.air.iii <- as.data.frame(data.Stat.3/(0.5*data$time.day[3]))
# Bind the data frames together row-wise
conc.air.i.2 <- rbind(conc.air.i, conc.air.ii)
# Average i and ii
conc.air.i.2 <- as.data.frame(t(colMeans(conc.air.i.2)))
# Bind the data frames together row-wise
conc.wb.air <- rbind(conc.air.i.2, conc.air.iii)
# Change column names of the last three columns
rownames(conc.wb.air) <- c("Conc.Air.1", "Conc.Air.2")
# Check total PCB
tPCB.conc.air <- rowSums(conc.wb.air, na.rm = TRUE)
# See
tPCB.conc.air

# (2) Wore WBs
# Read calculated average sampling rates
sr <- read.csv("Output/Data/csv/Ave.SRs.csv")
# Select only average sampling rate
sr <- sr[, 1:2]

# Select wore WBs ---------------------------------------------------------
{
  wb.Mi.o <- data[4, c(2, 3:175)]
  wb.Mi.h <- data[5, c(2, 3:175)]
  wb.Ea.h <- data[6, c(2, 3:175)]
  wb.Ea.o <- data[7, c(2, 3:175)]
  wb.Ya.o <- data[8, c(2, 3:175)]
  wb.Ya.h <- data[9, c(2, 3:175)]
  wb.An.h <- data[10, c(2, 3:175)]
  wb.An.o <- data[11, c(2, 3:175)]
  wb.Xu.h <- data[12, c(2, 3:175)]
  wb.Xu.o <- data[13, c(2, 3:175)]
  # Combined wore WBs
  wb.wr <- rbind(wb.Mi.o, wb.Mi.h, wb.Ea.o, wb.Ea.h, wb.Ya.o, wb.Ya.h,
                 wb.An.o, wb.An.h, wb.Xu.o, wb.Xu.h)
}

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
# Alternative method -> use a constant sampling rate
# wb_div_sr <- sweep(wb_common, 2, 1.5, FUN = "/")
# Extract time.day from wb.wr
wb_time_day <- wb.wr$time.day
# Divide wb_div_sr further by corresponding value in wb_time_day
conc.wb.wr <- sweep(wb_div_sr, 1, wb_time_day, FUN = "/")
rownames(conc.wb.wr) <- c('wb.Mi.o', 'wb.Mi.h', 'wb.Ea.o', 'wb.Ea.h',
                          'wb.Ya.o', 'wb.Ya.h', 'wb.An.o', 'wb.An.h',
                          'wb.Xu.o', 'wb.Xu.h')
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
prof.wb.air.conc[, 2:3] <- lapply(prof.wb.air.conc[, 2:3], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.air.conc$congener <- factor(prof.wb.air.conc$congener,
                                levels = unique(prof.wb.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.air.conc[, 2:3], na.rm = TRUE)

# (2) Wore WBs
tmp.wb.wr <- rowSums(conc.wb.wr, na.rm = TRUE)
prof.wb.wr.conc <- sweep(conc.wb.wr, 1, tmp.wb.wr, FUN = "/")
prof.wb.wr.conc <- data.frame(t(prof.wb.wr.conc))
congener <- rownames(prof.wb.wr.conc)
prof.wb.wr.conc <- cbind(congener, prof.wb.wr.conc)
rownames(prof.wb.wr.conc) <- NULL
prof.wb.wr.conc[, 2:11] <- lapply(prof.wb.wr.conc[, 2:11], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.wr.conc$congener <- factor(prof.wb.wr.conc$congener,
                               levels = unique(prof.wb.wr.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.wr.conc[, 2:11], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# Mi (1)
prof_combined.Mi <- prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Mi.o) %>%
              rename(Conc = wb.Mi.o) %>%
              mutate(Source = "Vol. 1 of")) %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Mi.h) %>%
              rename(Conc = wb.Mi.h) %>%
              mutate(Source = "Vol. 1 ho"))  # Change this label

# Plot
# Create the plot with the legend moved inside
p_prof_comb.Mi <- ggplot(prof_combined.Mi, aes(x = congener, y = Conc,
                                                 fill = Source)) +
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
                               "Vol. 1 of" = "#009E73",
                               "Vol. 1 ho" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Mi)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol1.png",
       plot = p_prof_comb.Mi, width = 10, height = 5, dpi = 500)

# Ea
prof_combined.Ea <- prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ea.o) %>%
              rename(Conc = wb.Ea.o) %>%
              mutate(Source = "Vol. 2 of")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ea.h) %>%
              rename(Conc = wb.Ea.h) %>%
              mutate(Source = "Vol. 2 ho"))

# Plots
p_prof_comb.Ea <- ggplot(prof_combined.Ea, aes(x = congener, y = Conc,
                                                 fill = Source)) +
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
                               "Vol. 2 of" = "#009E73",
                               "Vol. 2 ho" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Ea)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol2.png",
       plot = p_prof_comb.Ea, width = 10, height = 5, dpi = 500)

# Ya
prof_combined.Ya <- prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ya.o) %>%
              rename(Conc = wb.Ya.o) %>%
              mutate(Source = "Vol. 3 of")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ya.h) %>%
              rename(Conc = wb.Ya.h) %>%
              mutate(Source = "Vol. 3 ho"))

# Plots
p_prof_comb.Ya <- ggplot(prof_combined.Ya, aes(x = congener, y = Conc, fill = Source)) +
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
                               "Vol. 3 of" = "#009E73",
                               "Vol. 3 ho" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Ya)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol3.png",
       plot = p_prof_comb.Ya, width = 10, height = 5, dpi = 500)

# An
prof_combined.An <- prof.wb.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.An.o) %>%
              rename(Conc = wb.An.o) %>%
              mutate(Source = "Vol. 4 of")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.An.h) %>%
              rename(Conc = wb.An.h) %>%
              mutate(Source = "Vol. 4 ho"))

p_prof_comb.An <- ggplot(prof_combined.An, aes(x = congener, y = Conc,
                                                 fill = Source)) +
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
                               "Vol. 4 of" = "#009E73",
                               "Vol. 4 ho" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.An)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol4.png",
       plot = p_prof_comb.An, width = 10, height = 5, dpi = 500)

# Xu
prof_combined.Xu <-  prof.wb.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Xu.o) %>%
              rename(Conc = wb.Xu.o) %>%
              mutate(Source = "Vol. 5 of")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Xu.o) %>%
              rename(Conc = wb.Xu.o) %>%
              mutate(Source = "Vol. 5 ho"))

p_prof_comb.Xu <- ggplot(prof_combined.Xu, aes(x = congener, y = Conc,
                                                 fill = Source)) +
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
                               "Vol. 5 of" = "#009E73",
                               "Vol. 5 ho" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Xu)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol5.png",
       plot = p_prof_comb.Xu, width = 10, height = 5, dpi = 500)

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

# Ea
prof_combined_wide.Ea  <- prof_combined.Ea %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Ea[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Ea <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Ea

# Ya
prof_combined_wide.Ya  <- prof_combined.Ya %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Ya[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Ya <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Ya

# An
prof_combined_wide.An  <- prof_combined.An %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.An[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.An <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.An

# Xu
prof_combined_wide.Xu  <- prof_combined.Xu %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Xu[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Xu <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Xu

