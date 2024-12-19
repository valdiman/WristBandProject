# Concentration estimation

# Install packages
install.packages("readxl") #say no!
install.packages("ggplot2")

# Load libraries
{
  library(readxl)
  library(ggplot2)
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

# t test for tPCB ---------------------------------------------------------
# Checking normality
shapiro.test(log10(tPCB.conc.air)) # no needs but to be consistent
shapiro.test(log10(tPCB.conc.wb)) # needs to be log10 transformed

# two-sample t test
# (1) tPCB
tPCB.ttest <- t.test(log10(tPCB.conc.air), log10(tPCB.conc.wb))
print(tPCB.ttest) # No significant, p-value = 0.714

# Test if the difference is equal to 0
# Compute the differences
# Calculate differences for Amanda
diff_amanda.1 <- tPCB.conc.wb[1] - tPCB.conc.air[1]
diff_amanda.2 <- tPCB.conc.wb[2] - tPCB.conc.air[1]
# Calculate differences for Kay
diff_kay <- tPCB.conc.wb[3] - tPCB.conc.air[2]
# Calculate differences for Yau.weeks
diff_yau_1st <- tPCB.conc.wb[4] - tPCB.conc.air[3]
diff_yau_2nd <- tPCB.conc.wb[5] - tPCB.conc.air[4]
# Calculate differences for Yau.nw
diff_yau_nw <- tPCB.conc.wb[6] - tPCB.conc.air[5]
diff_yau_w <- tPCB.conc.wb[7] - tPCB.conc.air[5]
# Combine all differences into a single vector
all_differences <- c(diff_amanda.1, diff_amanda.2, diff_kay, diff_yau_1st,
                     diff_yau_2nd, diff_yau_nw, diff_yau_w)

# Print the combined differences
print(all_differences)
# Perform one-sample t-test
result <- t.test(all_differences, mu = 0)
# View the result
print(result) # Not significant, p-value = 0.975

# Using log10 transformation
log10.tPCB.conc.air <- log10(tPCB.conc.air)
log10.tPCB.conc.wb <- log10(tPCB.conc.wb)

# Test if the difference is equal to 0
# Compute the differences
# Calculate differences for Amanda
log10.diff_amanda.1 <- log10.tPCB.conc.wb[1] - log10.tPCB.conc.air[1]
log10.diff_amanda.2 <- log10.tPCB.conc.wb[2] - log10.tPCB.conc.air[1]
# Calculate differences for Kay
log10.diff_kay <- log10.tPCB.conc.wb[3] - log10.tPCB.conc.air[2]
# Calculate differences for Yau.weeks
log10.diff_yau_1st <- log10.tPCB.conc.wb[4] - log10.tPCB.conc.air[3]
log10.diff_yau_2nd <- log10.tPCB.conc.wb[5] - log10.tPCB.conc.air[4]
# Calculate differences for Yau.nw
log10.diff_yau_nw <- log10.tPCB.conc.wb[6] - log10.tPCB.conc.air[5]
log10.diff_yau_w <- log10.tPCB.conc.wb[7] - log10.tPCB.conc.air[5]
# Combine all differences into a single vector
log10.all_differences <- c(log10.diff_amanda.1, log10.diff_amanda.2,
                           log10.diff_kay, log10.diff_yau_1st,
                           log10.diff_yau_2nd, log10.diff_yau_nw,
                           log10.diff_yau_w)

# Print the combined differences
print(log10.all_differences)
# Perform one-sample t-test
# Not usre if this makes sense, using a mu = 0 when log10 is used
result <- t.test(log10.all_differences, mu = 0)
# View the result
print(result) # Not significant, p-value = 0.65

# t test for individual PCBs ----------------------------------------------
# Transpose the conc_air_common data frame
conc_air_common_transposed <- t(conc_air_common)
# Check normality
# Apply Shapiro-Wilk test to each column
shapiro_results.1 <- lapply(log10(conc.wb), shapiro.test)

# Extract p-values from Shapiro-Wilk test results
p_values.wb <- sapply(shapiro_results.1, function(x) x$p.value)

# Create data frame with congener names and p-values
p_values_wb <- data.frame(p_value.wb = p_values.wb)

# Count how many rows have p-values above 0.05
above_threshold_count <- sum(p_values.wb > 0.05, na.rm = TRUE)

# Print the count
print(above_threshold_count)

# Initialize an empty data frame with a single column for p-values and appropriate row names
shapiro_results.2 <- data.frame(p_value = rep(NA, ncol(conc_air_common_transposed)), 
                                row.names = colnames(conc_air_common_transposed))

# Populate the data frame with p-values
for (congener in colnames(conc_air_common_transposed)) {
  x <- log10(conc_air_common_transposed[, congener])
  if (sum(!is.na(x)) >= 3 && sum(is.na(x)) == 0) {
    result <- shapiro.test(x)
    shapiro_results.2[congener, "p_value"] <- result$p.value
  } else {
    shapiro_results.2[congener, "p_value"] <- NA
  }
}

# Replace 0 with NA in the p_value column
shapiro_results.2$p_value[shapiro_results.2$p_value == 0] <- NA

p_values_air <- shapiro_results.2

# Count how many rows have p-values above 0.05
above_threshold_count <- sum(p_values_air > 0.05, na.rm = TRUE)

# Print the count
print(above_threshold_count)

# Compute the differences. NA values are not included in this analysis.
# (1.1) Calculate differences for Amanda
diff_amanda.1 <- conc.wb[1, ] - conc_air_common_transposed[1, ]
diff_amanda.1 <- ifelse(is.na(diff_amanda.1) | is.na(conc_air_common_transposed[1, ]) | is.na(conc.wb[1, ]),
                        NA, diff_amanda.1)
# Change to data frame
diff_amanda.1 <- as.data.frame(diff_amanda.1)
# Add congener names to the columns
colnames(diff_amanda.1) <- colnames(conc_air_common_transposed)

# (1.2) Calculate differences for Amanda
diff_amanda.2 <- conc.wb[2, ] - conc_air_common_transposed[1, ]
diff_amanda.2 <- ifelse(is.na(diff_amanda.2) | is.na(conc_air_common_transposed[1, ]) | is.na(conc.wb[2, ]), NA, diff_amanda.2)
# Change to data frame
diff_amanda.2 <- as.data.frame(diff_amanda.2)
# Add congener names to the columns
colnames(diff_amanda.2) <- colnames(conc_air_common_transposed)

# (2) Calculate differences for Kay
diff_kay <- conc.wb[3, ] - conc_air_common_transposed[2, ]
diff_kay <- ifelse(is.na(diff_kay) | is.na(conc_air_common_transposed[2, ]) | is.na(conc.wb[3, ]), NA, diff_kay)
# Change to data frame
diff_kay <- as.data.frame(diff_kay)
# Add congener names to the columns
colnames(diff_kay) <- colnames(conc_air_common_transposed)

# (3.1) Calculate differences for Yau.weeks
diff_yau_1st <- conc.wb[4, ] - conc_air_common_transposed[3, ]
diff_yau_1st <- ifelse(is.na(diff_yau_1st) | is.na(conc_air_common_transposed[3, ]) | is.na(conc.wb[4, ]), NA, diff_yau_1st)
# Change to data frame
diff_yau_1st <- as.data.frame(diff_yau_1st)
# Add congener names to the columns
colnames(diff_yau_1st) <- colnames(conc_air_common_transposed)

# (3.2) Calculate differences for Yau.weeks
diff_yau_2nd <- conc.wb[5, ] - conc_air_common_transposed[4, ]
diff_yau_2nd <- ifelse(is.na(diff_yau_2nd) | is.na(conc_air_common_transposed[4, ]) | is.na(conc.wb[5, ]), NA, diff_yau_2nd)
# Change to data frame
diff_yau_2nd <- as.data.frame(diff_yau_2nd)
# Add congener names to the columns
colnames(diff_yau_2nd) <- colnames(conc_air_common_transposed)

# (4.1) Calculate differences for Yau.nw
diff_yau_nw <- conc.wb[6, ] - conc_air_common_transposed[5, ]
diff_yau_nw <- ifelse(is.na(diff_yau_nw) | is.na(conc_air_common_transposed[5, ]) | is.na(conc.wb[6, ]), NA, diff_yau_nw)
# Change to data frame
diff_yau_nw <- as.data.frame(diff_yau_nw)
# Add congener names to the columns
colnames(diff_yau_nw) <- colnames(conc_air_common_transposed)

# (4.2) Calculate differences for Yau.w
diff_yau_w <- conc.wb[7, ] - conc_air_common_transposed[5, ]
diff_yau_w <- ifelse(is.na(diff_yau_w) | is.na(conc_air_common_transposed[5, ]) | is.na(conc.wb[7, ]), NA, diff_yau_w)
# Change to data frame
diff_yau_w <- as.data.frame(diff_yau_w)
# Add congener names to the columns
colnames(diff_yau_w) <- colnames(conc_air_common_transposed)

# Combine all differences
all_differences <- rbind(diff_amanda.1, diff_amanda.2, diff_kay, diff_yau_1st,
                     diff_yau_2nd, diff_yau_nw, diff_yau_w)

# Create an empty data frame to store p-values
p_values <- data.frame(PCB = colnames(all_differences), p_value = NA)

# Iterate over each column of all_differences
for (i in seq_along(colnames(all_differences))) {
  col <- colnames(all_differences)[i]
  # Perform t-test for the current column
  t_test_result <- try(t.test(all_differences[[col]], mu = 0), silent = TRUE)
  
  # Extract the p-value if t-test was successful
  if (!inherits(t_test_result, "try-error")) {
    p_value <- t_test_result$p.value
  } else {
    # If t-test was not performed due to insufficient data or error, store NA
    p_value <- NA
  }
  
  # Store the p-value in the data frame
  p_values$p_value[i] <- p_value
}

# Plot p-values of PCBs ---------------------------------------------------
# Convert data frame to long format
data_long <- reshape2::melt(p_values)

# Reorder levels of the "PCB" column
data_long$PCB <- factor(data_long$PCB, levels = unique((p_values$PCB)))

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
p <- ggplot(data_long, aes(x = PCB, y = -log10(value))) +
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


