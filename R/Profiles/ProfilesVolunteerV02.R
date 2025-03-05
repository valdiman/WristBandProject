## Script to fabricate PCB profiles &
# and calculate cosine theta for similarity analysis

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
sr <- read.csv("Output/Data/csv/SamplingRates/Personal/PersonalAveSRV01.csv")
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
# Check sum of all PCBs (i.e. = 1)
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
# Check sum of all PCBs (i.e. = 1)
colSums(prof.wb.wr.conc[, 2:11], na.rm = TRUE)

# Concentration profile plots ---------------------------------------------
# Mi (1)
prof_combined.Mi <- prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Mi.o) %>%
              rename(Conc = wb.Mi.o) %>%
              mutate(Source = "Vol. 1 office-only")) %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Mi.h) %>%
              rename(Conc = wb.Mi.h) %>%
              mutate(Source = "Vol. 1 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 1 office-only",
                                            "Vol. 1 full-day")))

# Plot
# Create the plot with the legend moved inside
p_prof_comb.Mi <- ggplot(prof_combined.Mi, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +
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
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 1 office-only" = "#009E73",
                               "Vol. 1 full-day" = "#E69F00"),
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
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ea.o) %>%
              rename(Conc = wb.Ea.o) %>%
              mutate(Source = "Vol. 3 office-only")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ea.h) %>%
              rename(Conc = wb.Ea.h) %>%
              mutate(Source = "Vol. 3 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 3 office-only",
                                            "Vol. 3 full-day")))

# Plots
p_prof_comb.Ea <- ggplot(prof_combined.Ea, aes(x = congener, y = Conc,
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
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 3 office-only" = "#009E73",
                               "Vol. 3 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Ea)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol3.png",
       plot = p_prof_comb.Ea, width = 10, height = 5, dpi = 500)

# Ya
prof_combined.Ya <- prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ya.o) %>%
              rename(Conc = wb.Ya.o) %>%
              mutate(Source = "Vol. 2 office-only")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ya.h) %>%
              rename(Conc = wb.Ya.h) %>%
              mutate(Source = "Vol. 2 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 2 office-only",
                                            "Vol. 2 full-day")))

# Plots
p_prof_comb.Ya <- ggplot(prof_combined.Ya, aes(x = congener, y = Conc, fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",
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
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 2 office-only" = "#009E73",
                               "Vol. 2 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Ya)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol2.png",
       plot = p_prof_comb.Ya, width = 10, height = 5, dpi = 500)

# An
prof_combined.An <- prof.wb.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.An.o) %>%
              rename(Conc = wb.An.o) %>%
              mutate(Source = "Vol. 9 office-only")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.An.h) %>%
              rename(Conc = wb.An.h) %>%
              mutate(Source = "Vol. 9 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 9 office-only",
                                            "Vol. 9 full-day")))

p_prof_comb.An <- ggplot(prof_combined.An, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +
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
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 9 office-only" = "#009E73",
                               "Vol. 9 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.An)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol9.png",
       plot = p_prof_comb.An, width = 10, height = 5, dpi = 500)

# Xu
prof_combined.Xu <-  prof.wb.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Xu.o) %>%
              rename(Conc = wb.Xu.o) %>%
              mutate(Source = "Vol. 8 office-only")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Xu.h) %>%
              rename(Conc = wb.Xu.h) %>%
              mutate(Source = "Vol. 8 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 8 office-only",
                                            "Vol. 8 full-day")))

p_prof_comb.Xu <- ggplot(prof_combined.Xu, aes(x = congener, y = Conc,
                                                 fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +
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
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 8 office-only" = "#009E73",
                               "Vol. 8 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.Xu)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Vol8.png",
       plot = p_prof_comb.Xu, width = 10, height = 5, dpi = 500)


# Volunteers in location 1 (Conc.Air.1) Vol. 2 out
prof_combined.1 <-  prof.wb.air.conc %>%
  select(congener, Conc.Air.1) %>%
  rename(Conc = Conc.Air.1) %>%
  mutate(Source = "Air PCB office 1") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Ea.h) %>%
              rename(Conc = wb.Ea.h) %>%
              mutate(Source = "Vol. 3 full-day")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Mi.h) %>%
              rename(Conc = wb.Mi.h) %>%
              mutate(Source = "Vol. 1 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 1", "Vol. 1 full-day",
                                            "Vol. 3 full-day")))

p_prof_comb.1 <- ggplot(prof_combined.1, aes(x = congener, y = Conc,
                                               fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +
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
  scale_fill_manual(values = c("Air PCB office 1" = "blue",
                               "Vol. 1 full-day" = "#009E73",
                               "Vol. 3 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.92, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.1)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Office1.png",
       plot = p_prof_comb.1, width = 10, height = 5, dpi = 500)

# Volunteers in location 2 (Conc.Air.2)
prof_combined.2 <-  prof.wb.air.conc %>%
  select(congener, Conc.Air.2) %>%
  rename(Conc = Conc.Air.2) %>%
  mutate(Source = "Air PCB office 2") %>%  # Change this label
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.An.h) %>%
              rename(Conc = wb.An.h) %>%
              mutate(Source = "Vol. 9 full-day")) %>%
  bind_rows(prof.wb.wr.conc %>%
              select(congener, wb.Xu.h) %>%
              rename(Conc = wb.Xu.h) %>%
              mutate(Source = "Vol. 8 full-day")) %>%
  mutate(Source = factor(Source, levels = c("Air PCB office 2", "Vol. 8 full-day",
                                            "Vol. 9 full-day")))

p_prof_comb.2 <- ggplot(prof_combined.2, aes(x = congener, y = Conc,
                                             fill = Source)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 1, 
           color = "black",  # Add black edges to the bars
           linewidth = 0.2) +
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
  scale_fill_manual(values = c("Air PCB office 2" = "blue",
                               "Vol. 8 full-day" = "#009E73",
                               "Vol. 9 full-day" = "#E69F00"),
                    guide = guide_legend(key.size = unit(0.5, "lines"))) +  # Smaller legend squares
  theme(legend.position = c(0.93, 0.8),  # Inside the plot
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),  # Removes the legend title
        legend.text = element_text(size = 8, face = "bold"))

# Print the plots
print(p_prof_comb.2)

# Save plot in folder
ggsave("Output/Plots/Profiles/OfficeHome/prof_combined.Office2.png",
       plot = p_prof_comb.2, width = 10, height = 5, dpi = 500)

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



