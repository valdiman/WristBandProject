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
{
  data <- data.frame(read_excel("Data/Volunteers.xlsx", sheet = "Sheet1",
                                col_names = TRUE, col_types = NULL))
  data.2 <- data.frame(read_excel("Data/VolunteersV02.xlsx", sheet = "Sheet1",
                                  col_names = TRUE, col_types = NULL))
}

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
  # For Gift
  mass.Gi.Stat <- data.frame(t(colMeans(data[22:23, 4:176])))
  # For Xue
  mass.Xue.Stat <- data.2[3, 3:175]
}

# Mass accumulated in wore WBs --------------------------------------------
{
  wb.Mi.l <- data[1, c(2, 4:176)]
  wb.Mi.r <- data[2, c(2, 4:176)]
  wb.Ya.l <- data[3, c(2, 4:176)]
  wb.Ya.r <- data[4, c(2, 4:176)]
  wb.Ea.l <- data[5, c(2, 4:176)]
  wb.Ea.r <- data[6, c(2, 4:176)]
  wb.Cr.l <- data[9, c(2, 4:176)]
  wb.Cr.r <- data[10, c(2, 4:176)]
  wb.Hu.l <- data[11, c(2, 4:176)]
  wb.Hu.r <- data[12, c(2, 4:176)]
  wb.Xu.l <- data[17, c(2, 4:176)]
  wb.Xu.r <- data[18, c(2, 4:176)]
  wb.Gi.l <- data[25, c(2, 4:176)]
  wb.Gi.r <- data[24, c(2, 4:176)]
  wb.Xue.l <- data.2[13, 2:175]
  # Combined wore WBs
  wb.wr <- rbind(wb.Mi.l, wb.Mi.r, wb.Ya.l, wb.Ya.r, wb.Ea.l, wb.Ea.r,
                 wb.Cr.l, wb.Cr.r, wb.Hu.l, wb.Hu.r, wb.Xu.l, wb.Xu.r,
                 wb.Gi.l, wb.Gi.r, wb.Xue.l)
  # Add row names
  rownames(wb.wr) <- c("Mass.Mi.l.wr", "Mass.Mi.r.wr", "Mass.Ya.l.wr", "Mass.Ya.r.wr",
                       "Mass.Ea.l.wr", "Mass.Ea.r.wr", "Mass.Cr.l.wr", "Mass.Cr.r.wr",
                       "Mass.Hu.l.wr", "Mass.Hu.r.wr", "Mass.Xu.l.wr", "Mass.Xu.r.wr",
                       "Mass.Gi.l.wr", "Mass.Gi.r.wr", "Mass.Xue.l.wr")
}

# Estimate air concentrations ---------------------------------------------
# (1) Air WBs
# = massWB/(0.5*time.day)
conc.Mi.Ya <- as.data.frame(t(mass.Mi.Ya.Stat/(0.5*data$time.day[8])))
conc.Ea <- as.data.frame(t(mass.Ea.Stat/(0.5*data$time.day[7])))
conc.Cr <- as.data.frame(t(mass.Cr.Stat/(0.5*data$time.day[14])))
conc.Hu <- as.data.frame(t(mass.Hu.Stat/(0.5*data$time.day[13])))
conc.Xu <- as.data.frame(t(mass.Xu.Stat/(0.5*data$time.day[19])))
conc.Gi <- as.data.frame(t(mass.Gi.Stat/(0.5*data$time.day[22])))
conc.Xue <- as.data.frame(t(mass.Xue.Stat/(0.5*data.2$time.day[3])))

# Combine concentrations
conc.air <- cbind(conc.Mi.Ya, conc.Ea, conc.Cr, conc.Hu, conc.Xu, conc.Gi,
                  conc.Xue)
# Change column names of the last three columns
colnames(conc.air) <- c("Conc.Air.Mi.Ya", "Conc.Air.Ea", "Conc.Air.Cr",
                        "Conc.Air.Hu", "Conc.Air.Xu", "Conc.Air.Gi",
                        "Conc.Air.Xue")

# Check total PCB
tPCB.conc.air <- colSums(conc.air, na.rm = TRUE)
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
wb_common <- wb.wr[, intersect(colnames(wb.wr[-1]), congener_names)]
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
conc.wb <- sweep(wb_div_sr, 1, wb_time_day, FUN = "/")
rownames(conc.wb) <- c('wb.Mi.l', 'wb.Mi.r', 'wb.Ya.l', 'wb.Ya.r', 'wb.Ea.l',
                       'wb.Ea.r', 'wb.Cr.l', 'wb.Cr.r', 'wb.Hu.l', 'wb.Hu.r',
                       'wb.Xu.l', 'wb.Xu.r', 'wb.Gi.l', 'wb.Gi.r', 'wb.Xue.l')
# Check total PCB
tPCB.conc.wb <- rowSums(conc.wb, na.rm = TRUE)
# See
tPCB.conc.wb

# Transpose conc.wb
conc.air.t <- data.frame(t(conc.air))
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb), rownames(conc.air))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.air))
# Subset conc.air to include only the rows with matching row names
conc.air <- conc.air[matching_indices, ]
# Transpose conc.wb.air.t again
conc.air.t <- data.frame(t(conc.air))

# Profiles
# (1) Air WBs
tmp.wb.air <- colSums(conc.air, na.rm = TRUE)
prof.air.conc <- sweep(conc.air, 2, tmp.wb.air, FUN = "/")
#prof.wb.air.conc <- data.frame(t(prof.wb.air.conc))
congener <- rownames(prof.air.conc)
prof.air.conc <- cbind(congener, prof.air.conc)
rownames(prof.air.conc) <- NULL
prof.air.conc[, 2:8] <- lapply(prof.air.conc[, 2:8], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.air.conc$congener <- factor(prof.air.conc$congener,
                                levels = unique(prof.air.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.air.conc[, 2:8], na.rm = TRUE)

# (2) Wore WBs
tmp.wb.wr <- rowSums(conc.wb, na.rm = TRUE)
prof.wb.conc <- sweep(conc.wb, 1, tmp.wb.wr, FUN = "/")
prof.wb.conc <- data.frame(t(prof.wb.conc))
congener <- rownames(prof.wb.conc)
prof.wb.conc <- cbind(congener, prof.wb.conc)
rownames(prof.wb.conc) <- NULL
prof.wb.conc[, 2:16] <- lapply(prof.wb.conc[, 2:16], as.numeric)
#Then turn it back into a factor with the levels in the correct order
prof.wb.conc$congener <- factor(prof.wb.conc$congener,
                               levels = unique(prof.wb.conc$congener))
# Check sum of all PCBs (i.e., = 1)
colSums(prof.wb.conc[, 2:16], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# Mi (1)
prof_combined.Mi <- prof.air.conc %>%
  select(congener, Conc.Air.Mi.Ya) %>%
  rename(Conc = Conc.Air.Mi.Ya) %>%
  mutate(Source = "Air PCB") %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Mi.l) %>%
              rename(Conc = wb.Mi.l) %>%
              mutate(Source = "Vol. 1 nd")) %>%  # Change this label
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Mi.r) %>%
              rename(Conc = wb.Mi.r) %>%
              mutate(Source = "Vol. 1 d"))  # Change this label

# Plot
# Create the plot with the legend moved inside
p_prof_comb.Mi <- ggplot(prof_combined.Mi, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
                               "Vol. 1 nd" = "#009E73",
                               "Vol. 1 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Mi)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol1.png", plot = p_prof_comb.Mi,
       width = 10, height = 5, dpi = 500)

# Ya
prof_combined.Ya <- prof.air.conc %>%
  select(congener, Conc.Air.Mi.Ya) %>%
  rename(Conc = Conc.Air.Mi.Ya) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Ya.l) %>%
              rename(Conc = wb.Ya.l) %>%
              mutate(Source = "Vol. 2 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Ya.r) %>%
              rename(Conc = wb.Ya.r) %>%
              mutate(Source = "Vol. 2 nd"))

# Plots
p_prof_comb.Ya <- ggplot(prof_combined.Ya, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
print(p_prof_comb.Ya)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol2.png", plot = p_prof_comb.Ya,
       width = 10, height = 5, dpi = 500)

# Ea
prof_combined.Ea <- prof.air.conc %>%
  select(congener, Conc.Air.Ea) %>%
  rename(Conc = Conc.Air.Ea) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Ea.l) %>%
              rename(Conc = wb.Ea.l) %>%
              mutate(Source = "Vol. 3 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Ea.r) %>%
              rename(Conc = wb.Ea.r) %>%
              mutate(Source = "Vol. 3 d"))

# Plots
p_prof_comb.Ea <- ggplot(prof_combined.Ea, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the barshttp://127.0.0.1:8373/graphics/plot_zoom_png?width=1872&height=900
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
print(p_prof_comb.Ea)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol3.png", plot = p_prof_comb.Ea,
       width = 10, height = 5, dpi = 500)

# Cr
prof_combined.Cr <- prof.air.conc %>%
  select(congener, Conc.Air.Cr) %>%
  rename(Conc = Conc.Air.Cr) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Cr.l) %>%
              rename(Conc = wb.Cr.l) %>%
              mutate(Source = "Vol. 4 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Cr.r) %>%
              rename(Conc = wb.Cr.r) %>%
              mutate(Source = "Vol. 4 d"))

p_prof_comb.Cr <- ggplot(prof_combined.Cr, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)http://127.0.0.1:8373/graphics/plot_zoom_png?width=1856&height=861
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
print(p_prof_comb.Cr)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol4.png", plot = p_prof_comb.Cr,
       width = 10, height = 5, dpi = 500)

# Hu
prof_combined.Hu <- prof.air.conc %>%
  select(congener, Conc.Air.Hu) %>%
  rename(Conc = Conc.Air.Hu) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Hu.l) %>%
              rename(Conc = wb.Hu.l) %>%
              mutate(Source = "Vol. 5 d")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Hu.r) %>%
              rename(Conc = wb.Hu.r) %>%
              mutate(Source = "Vol. 5 nd"))

p_prof_comb.Hu <- ggplot(prof_combined.Hu, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
print(p_prof_comb.Hu)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol5.png", plot = p_prof_comb.Hu,
       width = 10, height = 5, dpi = 500)

# Xu
prof_combined.Xu <- prof.air.conc %>%
  select(congener, Conc.Air.Xu) %>%
  rename(Conc = Conc.Air.Xu) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Xu.l) %>%
              rename(Conc = wb.Xu.l) %>%
              mutate(Source = "Vol. 6 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Xu.r) %>%
              rename(Conc = wb.Xu.r) %>%
              mutate(Source = "Vol. 6 d"))

# Plots
p_prof_comb.Xu <- ggplot(prof_combined.Xu, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
print(p_prof_comb.Xu)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol6.png", plot = p_prof_comb.Xu,
       width = 10, height = 5, dpi = 500)

# Gift
prof_combined.Gi <- prof.air.conc %>%
  select(congener, Conc.Air.Gi) %>%
  rename(Conc = Conc.Air.Gi) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Gi.l) %>%
              rename(Conc = wb.Gi.l) %>%
              mutate(Source = "Vol. 7 nd")) %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Gi.r) %>%
              rename(Conc = wb.Gi.r) %>%
              mutate(Source = "Vol. 7 d"))

# Plots
p_prof_comb.Gi <- ggplot(prof_combined.Gi, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
  xlab("") +
  ylim(0, 0.7) +
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
                               "Vol. 7 nd" = "#009E73",
                               "Vol. 7 d" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Gi)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol7.png", plot = p_prof_comb.Gi,
       width = 10, height = 5, dpi = 500)

# Xuefang
prof_combined.Xue <- prof.air.conc %>%
  select(congener, Conc.Air.Xue) %>%
  rename(Conc = Conc.Air.Xue) %>%
  mutate(Source = "Air PCB") %>%
  bind_rows(prof.wb.conc %>%
              select(congener, wb.Xue.l) %>%
              rename(Conc = wb.Xue.l) %>%
              mutate(Source = "Vol. 8 nd"))

# Plots
p_prof_comb.Xue <- ggplot(prof_combined.Xue, aes(x = congener, y = Conc,
                                                  fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +  # Set the thickness of the black edges (fine line)
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
                               "Vol. 8 nd" = "#009E73"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Xue)

# Save plot in folder
ggsave("Output/Plots/Profiles/prof_combined.Vol8.png", plot = p_prof_comb.Xue,
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

# Gift
prof_combined_wide.Gi  <- prof_combined.Gi %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Gi[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Gi <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Gi

# Xuefang
prof_combined_wide.Xue  <- prof_combined.Xue %>%
  pivot_wider(names_from = Source, values_from = Conc)
# Create a matrix from the concentration values (excluding the congener column)
cosine_matrix <- as.matrix(prof_combined_wide.Xue[, -1])
# Calculate cosine similarity between the columns (now rows after transposing)
cosine_similarity.Xue <- cosine(cosine_matrix)
# View the resulting cosine similarity matrix
cosine_similarity.Xue

