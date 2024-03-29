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
# Use cte SR (0.5 and 1)
# wb_div_sr <- wb_common/1
# Divide further by corresponding value in wb.5d$time.day
conc.wb <- wb_div_sr/wb_time_day
rownames(conc.wb) <- c('wb.amanda.r', 'wb.amanda.l', 'wb.kay', 'wb.yau.1st',
                       'wb.yau.2nd', 'wb.yau.nw', 'wb.yau.w')

# Match congeners in both dataset -----------------------------------------
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb), rownames(conc.air))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.air))
# Subset conc.air to include only the rows with matching row names
conc_air_common <- conc.air[matching_indices, ]

# Sum PCBs ----------------------------------------------------------------
tPCB.conc.wb <- rowSums(conc.wb, na.rm = TRUE)
# Convert list column to numeric, replacing NA values with 0
conc_air_common_numeric <- as.data.frame(sapply(conc_air_common, as.numeric))
tPCB.conc.air <- colSums(conc_air_common_numeric, na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Combine data for error calculation and plotting -------------------------
# Create data frames
{
  # For Amanda
  df_amanda_r <- data.frame(
    Conc.Air.Amanda = tPCB.conc.air["Conc.Air.Amanda"],
    wb.amanda.r = tPCB.conc.wb["wb.amanda.r"]
  )
  
  df_amanda_l <- data.frame(
    Conc.Air.Amanda = tPCB.conc.air["Conc.Air.Amanda"],
    wb.amanda.l = tPCB.conc.wb["wb.amanda.l"]
  )
  
  # For Kay
  df_kay <- data.frame(
    Conc.Air.Kay = tPCB.conc.air["Conc.Air.Kay"],
    wb.amanda.r = tPCB.conc.wb["wb.kay"]
  )
  
  # For Yau
  df_yau_1st <- data.frame(
    Conc.Air.Yau.1st = tPCB.conc.air["Conc.Air.Yau.1st"],
    wb.yau.1st = tPCB.conc.wb["wb.yau.1st"]
  )
  
  df_yau_2nd <- data.frame(
    Conc.Air.Yau.2nd = tPCB.conc.air["Conc.Air.Yau.2nd"],
    wb.yau.2nd = tPCB.conc.wb["wb.yau.2nd"]
  )
  
  df_yau_nw <- data.frame(
    Conc.Air.Yau.w = tPCB.conc.air["Conc.Air.Yau.w"],
    wb.yau.nw = tPCB.conc.wb["wb.yau.nw"]
  )
  
  df_yau_w <- data.frame(
    Conc.Air.Yau.w = tPCB.conc.air["Conc.Air.Yau.w"],
    wb.yau.w = tPCB.conc.wb["wb.yau.w"]
  )
}

# Change column names
{
  colnames(df_amanda_r) <- c("Conc.Air", "Conc.WB")
  colnames(df_amanda_l) <- c("Conc.Air", "Conc.WB")
  colnames(df_kay) <- c("Conc.Air", "Conc.WB")
  colnames(df_yau_1st) <- c("Conc.Air", "Conc.WB")
  colnames(df_yau_2nd) <- c("Conc.Air", "Conc.WB")
  colnames(df_yau_nw) <- c("Conc.Air", "Conc.WB")
  colnames(df_yau_w) <- c("Conc.Air", "Conc.WB")
}

# Combine all data frames
combined_data <- rbind(df_amanda_r, df_amanda_l, df_kay, df_yau_1st, df_yau_2nd,
                       df_yau_nw, df_yau_w)

# Extract the group names from the row names
combined_data$Group <- gsub("^Conc\\.Air\\.|^wb\\.", "", rownames(combined_data))

# Change the group names
{
  combined_data$Group <- gsub("Amanda$", "A.r", combined_data$Group)
  combined_data$Group <- gsub("Amanda1$", "A.l", combined_data$Group)
  combined_data$Group <- gsub("Kay", "K", combined_data$Group)
  combined_data$Group <- gsub("Yau.1st$", "Y.1st", combined_data$Group)
  combined_data$Group <- gsub("Yau.2nd$", "Y.2nd", combined_data$Group)
  combined_data$Group <- gsub("Yau\\.w$", "Y.nw", combined_data$Group)
  combined_data$Group <- gsub("Yau\\.w1$", "Y.w", combined_data$Group)
}

# Estimate error (factor of 2) --------------------------------------------
# Estimate a factor of 2 between observations and predictions
combined_data$factor2 <- combined_data$Conc.WB/combined_data$Conc.Air

# Calculate the percentage of observations within the factor of 2
factor2_percentage <- nrow(combined_data[combined_data$factor2 > 0.5 & combined_data$factor2 < 2, ])/nrow(combined_data)*100

# Estimate RMSE and MAE ---------------------------------------------------
# log10 scale
rmse <- sqrt(mean((log10(combined_data$Conc.Air) - log10(combined_data$Conc.WB))^2))
mae <- mean(abs(log10(combined_data$Conc.Air) - log10(combined_data$Conc.WB)))
print(paste("RMSE:", rmse))
print(paste("MAE:", mae))

# Plot --------------------------------------------------------------------
# Choose a color palette from ColorBrewer
num_groups <- length(unique(combined_data$Group))
color_palette <- brewer.pal(num_groups, "Set1")

# Plot with different colors for each group
plotWB <- ggplot(combined_data, aes(x = Conc.Air, y = Conc.WB, color = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  scale_y_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Air Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  ylab(expression(bold("WB Concentration " *Sigma*"PCB (ng/m"^3*")"))) +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) +
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +
  scale_color_manual(values = color_palette) +
  labs(color = "Participants") 
# If text needs to be included in the plot
#annotate("text", x = 10000, y = 20, hjust = 1, vjust = 1,
#           label = expression("Sampling rate = 1.0 (m"^3*"/d)"),
#           size = 3, color = "black")

# Print the plot. Warning message. No worries
print(plotWB)

# Save plot in folder
ggsave("Output/Plots/tPCB.png", plot = plotWB, width = 6, height = 5, dpi = 500)

