# Concentration estimation

# Install packages
install.packages("readxl") #say no!
install.packages("ggplot2")
install.packages("lme4")
install.packages("Matrix")
install.packages("Matrix", type = "binary")
install.packages("lmerTest")
install.packages("tidyverse")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(lme4) # performs lme
  library(Matrix)
  library(lmerTest) # gets the p-value from lme
  library(tidyverse)
}

# Read measured values from excel -----------------------------------------
{
  data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                                       col_names = TRUE, col_types = NULL))
  data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                    col_names = TRUE, col_types = NULL))
  data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                    col_names = TRUE, col_types = NULL))
  data.yau2 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet2",
                                     col_names = TRUE, col_types = NULL)) 
}

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

# Normality check ---------------------------------------------------------
# Plots
# (1) Histograms
hist(tPCB.conc.air)
hist(log10(tPCB.conc.air))
hist(tPCB.conc.wb)
hist(log10(tPCB.conc.wb))
# (2) Normal Q-Q plot
qqnorm(tPCB.conc.air)
qqline(tPCB.conc.air)
qqnorm(log10(tPCB.conc.air))
qqline(log10(tPCB.conc.air))

qqnorm(tPCB.conc.wb)
qqline(tPCB.conc.wb)
qqnorm(log10(tPCB.conc.wb))
qqline(log10(tPCB.conc.wb))

# Shapiro test
shapiro.test(log10(tPCB.conc.air)) # no needs but to be consistent
shapiro.test(log10(tPCB.conc.wb)) # needs to be log10 transformed

# Merge data
# Create data frames for tPCB.conc.wb and tPCB.conc.air
tPCB.conc.wb_df <- data.frame(
  Description = c("1.wb.amanda.r", "1.wb.amanda.l", "2.wb.kay", "3.wb.yau.1st", "4.wb.yau.2nd", "5.wb.yau.nw", "5.wb.yau.w"),
  Value = c(606.2357, 289.1449, 104.8486, 186.6856, 242.9024, 187.1839, 166.6439)
)

# Add a new column "method" with all rows equal to 2
tPCB.conc.wb_df$method <- 2

tPCB.conc.air_df <- data.frame(
  Description = c("1.Conc.Air.Amanda", "2.Conc.Air.Kay", "3.Conc.Air.Yau.1st", "4.Conc.Air.Yau.2nd", "5.Conc.Air.Yau.w"),
  Value = c(289.0832, 351.7530, 160.5868, 230.2299, 224.2439)
)

# Add a new column "method" with all rows equal to 2
tPCB.conc.air_df$method <- 1

# Combine the datasets longitudinally
merged_df <- rbind(tPCB.conc.wb_df, tPCB.conc.air_df)

# Extract the first number from the left in the "Description" column
merged_df$ID <- as.integer(gsub("\\D*(\\d+).*", "\\1", merged_df$Description))

# Print the modified data frame
print(merged_df)

# Anova for tPCB ----------------------------------------------------------
# Anova (1)
fit = lmer(Value ~ method + (1|ID), data = merged_df)
summary(fit)
anova(fit) # p-value = 0.9981!

# Anova (2)
fit = lmer(log10(Value) ~ method + (1|ID), data = merged_df)
summary(fit)
anova(fit) # p-value = 0.6592!

# Anova for individual PCBs -----------------------------------------------
# Convert row names to a new column, typically named 'rowname'
conc_air_v02 <- tibble::rownames_to_column(conc_air_common, var = "ID")

# Reshape data
long_data_conc_air <- pivot_longer(
  conc_air_v02,
  cols = -ID,  # Select all columns except the ID for pivoting
  names_to = "Description",  # This will hold the original column names as descriptions
  values_to = "Value"  # This will hold the numerical values
)

# Add method 1
long_data_conc_air$method <- 1  # Add a method column with all values set to 1

unique_descriptions <- unique(long_data_conc_air$Description)

description_to_number <- data.frame(
  Description = unique_descriptions,
  Number = seq_along(unique_descriptions)
)

long_data_conc_air <- merge(long_data_conc_air, description_to_number,
                            by = "Description", all.x = TRUE)


long_data_conc_air$Description <- paste0(long_data_conc_air$Number, ".",
                                         long_data_conc_air$Description)

long_data_conc_air$Number <- NULL

conc_wb_v02 <- tibble::rownames_to_column(conc.wb, var = "Description")

long_data_conc_wb <- pivot_longer(
  conc_wb_v02,
  cols = -Description,    # Select all columns except the Description for pivoting
  names_to = "ID",        # This will hold the Congeners names (column names)
  values_to = "Value"     # This will hold the numerical values
)

long_data_conc_wb$method <- 2  # Add a method column with all values set to 2

# Manually add the front number based on the pattern you described
long_data_conc_wb$Description <- paste0(
  ifelse(grepl("wb\\.amanda", long_data_conc_wb$Description), "1", 
         ifelse(grepl("wb\\.kay", long_data_conc_wb$Description), "2",
                ifelse(grepl("wb\\.yau\\.1st", long_data_conc_wb$Description), "3",
                       ifelse(grepl("wb\\.yau\\.2nd", long_data_conc_wb$Description), "4",
                              ifelse(grepl("wb\\.yau\\.[nw]", long_data_conc_wb$Description), "5", NA)
                       )
                )
         )
  ),
  ".", long_data_conc_wb$Description
)

# If you have multiple occurrences of "wb.yau.nw" and "wb.yau.w", adjust them manually
long_data_conc_wb$Description[long_data_conc_wb$Description == "5.wb.yau.nw"] <- "5.wb.yau.nw"
long_data_conc_wb$Description[long_data_conc_wb$Description == "5.wb.yau.w"] <- "5.wb.yau.w"

# Change the column name from "ID" to "congener"
names(long_data_conc_wb)[names(long_data_conc_wb) == "ID"] <- "congener"
names(long_data_conc_air)[names(long_data_conc_air) == "ID"] <- "congener"

# Combine the datasets longitudinally
merged_df <- rbind(long_data_conc_wb, long_data_conc_air)

# Extract the first number from the left in the "Description" column
merged_df$ID <- as.integer(gsub("\\D*(\\d+).*", "\\1", merged_df$Description))

# Split the dataframe by congener
congener_list <- split(merged_df, merged_df$congener)

# Initialize a data frame to store the results
p_values <- data.frame(congener = character(), p_value = numeric())

# Loop over each congener dataframe
for (congener_df in congener_list) {
  # Fit the model
  fit <- lmer(Value ~ method + (1|ID), data = congener_df)
  
  # Get the ANOVA table
  anova_table <- anova(fit)
  
  # Extract the p-value
  p_value <- anova_table$'Pr(>F)'[1]
  
  # Store the congener name and p-value
  p_values <- rbind(p_values,
                    data.frame(congener = congener_df$congener[1],
                               p_value = p_value))
}

# Print the results
print(p_values)

# Calculate the percentage of p-values above 0.05
percentage_above_005 <- (sum(p_values$p_value > 0.05) / nrow(p_values)) * 100

# log10
# Initialize a data frame to store the results
p_values <- data.frame(congener = character(), p_value = numeric())

# Loop over each congener dataframe
for (congener_df in congener_list) {
  # Replace 0 and NA values with a small positive value
  congener_df$Value <- ifelse(congener_df$Value <= 0 | is.na(congener_df$Value),
                              1e-6, congener_df$Value)
  
  # Calculate log10 of the "Value" column
  congener_df$log_Value <- log10(congener_df$Value)
  
  # Fit the model
  fit <- lmer(log_Value ~ method + (1|ID), data = congener_df)
  
  # Get the ANOVA table
  anova_table <- anova(fit)
  
  # Extract the p-value
  p_value <- anova_table$'Pr(>F)'[1]
  
  # Store the congener name and p-value
  p_values <- rbind(p_values,
                    data.frame(congener = congener_df$congener[1],
                               p_value = p_value))
}

# Print the results
print(p_values)

# Calculate the percentage of p-values above 0.05
percentage_above_005 <- (sum(p_values$p_value > 0.05) / nrow(p_values)) * 100

# Plot p-values of PCBs ---------------------------------------------------
# Need to reorganize the PCB congenes to be plotted.

# Calculate the number above and below -log10(0.05)
# Calculate the number of values greater than -log10(0.05)
above_threshold_count <- sum(-log10(p_values$p_value) > -log10(0.05),
                             na.rm = TRUE)
below_threshold_count <- sum(-log10(p_values$p_value) <= -log10(0.05),
                             na.rm = TRUE)

# Calculate the total number of values
total_count <- length(na.omit(p_values$p_value))

# Calculate the fraction of values greater than -log10(0.05) relative to the total
above_threshold_percentage <- above_threshold_count / total_count
below_threshold_percentage <- below_threshold_count / total_count

# Calculate the number above and below the threshold and convert to percentage
above_threshold_percentage <- round((above_threshold_percentage * 100),
                                    digits = 0)
below_threshold_percentage <- round((below_threshold_percentage * 100),
                                    digits = 0)

# Create the ggplot
p <- ggplot(plot_data, aes(x = congener, y = -log10(value))) +
  geom_point(shape = 19, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  labs(x = NULL, y = "-log10(p-value)") +
  theme_bw() +
  theme(aspect.ratio = 6/18) +
  theme(axis.text.x = element_text(size = 6, angle = 60, hjust = 1)) +  # Rotate x-axis labels
  annotate("text", x = 145, y = 3.9,
           label = paste("% > threshold:", above_threshold_percentage),
           color = "blue", size = 4) +
  annotate("text", x = 145, y = 3.5,
           label = paste("% < threshold:", below_threshold_percentage),
           color = "blue", size = 4) +
  annotate("text", x = 62, y = -log10(0.035),
           label = paste("p-value = 0.05"), color = "red", size = 4)

# See plot
print(p)

# Save plot in folder
ggsave("Output/Plots/pvaluePCBi.png", plot = p, width = 10,
       height = 5, dpi = 500)


