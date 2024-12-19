# Concentration estimation

# Install packages
install.packages("readxl") #say no!
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("RColorBrewer")
install.packages("scales")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(scales)
}

# Read measured values from excel -----------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                  col_names = TRUE, col_types = NULL))
data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                  col_names = TRUE, col_types = NULL))
data.yau2 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet2",
                                   col_names = TRUE, col_types = NULL))

# Calculate air PCB concentration from static WBs -------------------------
# Select WB to estimate airborne concentration & WB wore 5 days
# (1) Amanda
{
  data.amanda.1 <- data.amanda[1:3, ]
  # Average 3 WBs
  data.amanda.2 <- colMeans(data.amanda.1[, 3:175])
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc.amamda <- as.data.frame(data.amanda.2/(0.5*data.amanda[1,1]))
  colnames(conc.amamda) <- "Conc.Air.Amanda"
}

# (2) Kay
{
  # Select WBs to calculate air concentration
  data.kay.1 <- data.kay[1:3, ]
  # Average 3 WBs
  data.kay.2 <- colMeans(data.kay.1[, 3:175])
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc.kay <- as.data.frame(data.kay.2/(0.5*data.kay[1,1]))
  colnames(conc.kay) <- "Conc.Air.Kay"
}

# (3) Ya'u
{
  # Select WBs to calculate air concentration
  data.yau.1st <- data.yau[3, 4:176]
  data.yau.2nd <- data.yau[6, 4:176]
  data.yau.w <- data.yau2[3, 5:177]
  # Calculate air concentration in ng/m3 for Ya'u
  conc.yau.1st <- as.data.frame(data.yau.1st / (0.5 * data.yau[3, 1]))
  conc.yau.2nd <- as.data.frame(data.yau.2nd / (0.5 * data.yau[6, 1]))
  conc.yau.w <- as.data.frame(data.yau.w / (0.5 * data.yau2[3, 1]))
  # Transpose data frames
  conc.yau.1st <- as.data.frame(t(conc.yau.1st))
  conc.yau.2nd <- as.data.frame(t(conc.yau.2nd))
  conc.yau.w <- as.data.frame(t(conc.yau.w))
}

# Combine concentrations
conc.air <- cbind(conc.amamda, conc.kay, conc.yau.1st, conc.yau.2nd, conc.yau.w)
# Change column names of the last three columns
colnames(conc.air)[(ncol(conc.air)-2):ncol(conc.air)] <- c("Conc.Air.Yau.1st",
                                                           "Conc.Air.Yau.2nd",
                                                           "Conc.Air.Yau.w")

# Read calculated average sampling rates ----------------------------------
sr <- read.csv("Output/Data/csv/Ave.SRs.csv")
# Select only average sampling rate
sr <- sr[, 1:2]

# Select wore WBs for 5 days and time -------------------------------------
{
  wb.amamda.5d.r <- data.amanda[8, c(1, 3:175)]
  wb.amamda.5d.l <- data.amanda[13, c(1, 3:175)]
  wb.kay.5d <- data.kay[8, c(1, 3:175)]
  wb.yau.5d.1st <- data.yau[9, c(1, 4:176)]
  wb.yau.5d.2nd <- data.yau[12, c(1, 4:176)]
  wb.yau.5d.nw <- data.yau2[6, c(1, 5:177)]
  wb.yau.5d.w <- data.yau2[9, c(1, 5:177)]
}
# Combined wore WB 5 days
wb.5d <- rbind(wb.amamda.5d.r, wb.amamda.5d.l, wb.kay.5d, wb.yau.5d.1st,
               wb.yau.5d.2nd, wb.yau.5d.nw, wb.yau.5d.w)

# Estimate air concentration using SR, mass of wore WBs & time ------------
# Extract congener names from sr
congener_names <- sr$congener
# Subset wb.5d to include only the common congeners
wb_common <- wb.5d[, intersect(colnames(wb.5d), congener_names)]
# Subset sr to include only the common congeners
sr_common <- sr[sr$congener %in% colnames(wb_common), ]
# Extract time.day from wb.5d
wb_time_day <- wb.5d$time.day
# Divide each element in wb_common by corresponding element in sr_common$Average_Sampling_Rate
wb_div_sr <- wb_common/sr_common$Average_Sampling_Rate
# Using cte sampling rate
# wb_div_sr <- wb_common/1
# Divide further by corresponding value in wb.5d$time.day
conc.wb <- wb_div_sr/wb_time_day
rownames(conc.wb) <- c('wb.amanda.r', 'wb.amanda.l', 'wb.kay', 'wb.yau.1st',
                       'wb.yau.2nd', 'wb.yau.nw', 'wb.yau.w')

# Transpose conc.wb.t
conc.wb.t <- as.data.frame(t(conc.wb))

# Match congeners in both dataset -----------------------------------------
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb), rownames(conc.air))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.air))
# Subset conc.air to include only the rows with matching row names
conc_air_common <- conc.air[matching_indices, ]

# Combine data for error calculation and plotting -------------------------
  # Create data frames
  {
    # For Amanda
    df_amanda_r <- data.frame(
      Conc.Air.Amanda = conc_air_common["Conc.Air.Amanda"],
      wb.amanda.r = conc.wb.t["wb.amanda.r"],
      Participant = "A.r"
    )
    
    df_amanda_l <- data.frame(
      Conc.Air.Amanda = conc_air_common["Conc.Air.Amanda"],
      wb.amanda.l = conc.wb.t["wb.amanda.l"],
      Participant = "A.l"
    )
    
    # For Kay
    df_kay <- data.frame(
      Conc.Air.Kay = conc_air_common["Conc.Air.Kay"],
      wb.amanda.r = conc.wb.t["wb.kay"],
      Participant = "K"
    )
    
    # For Yau
    df_yau_1st <- data.frame(
      Conc.Air.Yau.1st = conc_air_common["Conc.Air.Yau.1st"],
      wb.yau.1st = conc.wb.t["wb.yau.1st"],
      Participant = "Y.1st"
    )
    
    df_yau_2nd <- data.frame(
      Conc.Air.Yau.2nd = conc_air_common["Conc.Air.Yau.2nd"],
      wb.yau.2nd = conc.wb.t["wb.yau.2nd"],
      Participant = "Y.2nd"
    )
    
    df_yau_nw <- data.frame(
      Conc.Air.Yau.w = conc_air_common["Conc.Air.Yau.w"],
      wb.yau.nw = conc.wb.t["wb.yau.nw"],
      Participant = "Y.nw"
    )
    
    df_yau_w <- data.frame(
      Conc.Air.Yau.w = conc_air_common["Conc.Air.Yau.w"],
      wb.yau.w = conc.wb.t["wb.yau.w"],
      Participant = "Y.w"
    )
  }

# Change column names
{
  colnames(df_amanda_r) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_amanda_l) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_kay) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_yau_1st) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_yau_2nd) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_yau_nw) <- c("Conc.Air", "Conc.WB", "Participant")
  colnames(df_yau_w) <- c("Conc.Air", "Conc.WB", "Participant")
}

# Combine all data frames
combined_data <- rbind(df_amanda_r, df_amanda_l, df_kay, df_yau_1st, df_yau_2nd,
                       df_yau_nw, df_yau_w)

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
combined_data$factor2 <- combined_data$Conc.WB/combined_data$Conc.Air

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(combined_data[combined_data$factor2 > 0.5 & combined_data$factor2 < 2, ])/nrow(combined_data)*100

# Plot --------------------------------------------------------------------
# Filter out NA and 0 values
filtered_data <- combined_data[complete.cases(combined_data) & combined_data$Conc.Air > 0 & combined_data$Conc.WB > 0, ]

# Choose a color palette from ColorBrewer
num_groups <- length(unique(filtered_data$Participant))
color_palette <- brewer.pal(num_groups, "Set1")

# Plot with different colors for each group
plotWB <- ggplot(filtered_data, aes(x = Conc.Air, y = Conc.WB, color = Participant)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(aspect.ratio = 15/15,
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_y_log10(limits = c(0.00001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.00001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Air Concentration PCBi (ng/m"^3*")"))) +
  ylab(expression(bold("WB Concentration PCBi (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  scale_color_manual(values = color_palette) +
  labs(color = "Participants") +
  annotation_logticks(sides = "bl") +
  annotate("text", x = 10000, y = 0.00001, hjust = 1, vjust = 1,
           label = expression("Sampling rate = 1.0 (m"^3*"/d)"),
           size = 5, color = "black")

# Print the plots
print(plotWB)

# Save plot in folder
ggsave("Output/Plots/PCBiSR=1.png", plot = plotWB, width = 10, height = 10, dpi = 500)
