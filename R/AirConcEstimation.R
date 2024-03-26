# Concentration estimation

# Install packages
install.packages("readxl") #say no!
install.packages("ggplot2")
install.packages("gridExtra")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
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
data.amanda.1 <- data.amanda[1:3, ]
# Average 3 WBs
data.amanda.2 <- colMeans(data.amanda.1[, 3:175])
# Calculate air concentration in ng/m3
# = massWB/(0.5*time.day)
conc.amamda <- as.data.frame(data.amanda.2/(0.5*data.amanda[1,1]))
colnames(conc.amamda) <- "Conc.Air.Amanda"

# (2) Kay
# Select WBs to calculate air concentration
data.kay.1 <- data.kay[1:3, ]
# Average 3 WBs
data.kay.2 <- colMeans(data.kay.1[, 3:175])
# Calculate air concentration in ng/m3
# = massWB/(0.5*time.day)
conc.kay <- as.data.frame(data.kay.2/(0.5*data.kay[1,1]))
colnames(conc.kay) <- "Conc.Air.Kay"

# (3) Ya'u
# Select WBs to calculate air concentration
data.yau.1st <- data.yau[3, 4:176]
data.yau.2nd <- data.yau[6, 4:176]
data.yau.w <- data.yau2[3, 5:177]

# Calculate air concentration in ng/m3
conc.yau.1st <- as.data.frame(matrix(data.yau.1st/(0.5*data.yau[3, 1]), ncol = 1))
rownames(conc.yau.1st) <- colnames(data.yau.1st)
colnames(conc.yau.1st) <- "Conc.Air.Yau.1st"
conc.yau.2nd <- as.data.frame(matrix(data.yau.2nd/(0.5*data.yau[6, 1]), ncol = 1))
rownames(conc.yau.2nd) <- colnames(data.yau.2nd)
colnames(conc.yau.2nd) <- "Conc.Air.Yau.2nd"
conc.yau.w <- as.data.frame(matrix(data.yau.w/(0.5*data.yau2[3, 1]), ncol = 1))
rownames(conc.yau.w) <- colnames(data.yau.w)
colnames(conc.yau.w) <- "Conc.Air.Yau.w"
# Combine concentrations
conc.air <- cbind(conc.amamda, conc.kay, conc.yau.1st, conc.yau.2nd, conc.yau.w)

# Read calculated average sampling rates ----------------------------------
sr <- read.csv("Output/Data/csv/Ave.SRs.csv")
# Select only average sampling rate
sr <- sr[, 1:2]

# Select wore WBs for 5 days and time -------------------------------------
wb.amamda.5d.r <- data.amanda[8, c(1, 3:175)]
wb.amamda.5d.l <- data.amanda[13, c(1, 3:175)]
wb.kay.5d <- data.kay[8, c(1, 3:175)]
wb.yau.5d.1st <- data.yau[9, c(1, 4:176)]
wb.yau.5d.2nd <- data.yau[12, c(1, 4:176)]
wb.yau.5d.nw <- data.yau2[6, c(1, 5:177)]
wb.yau.5d.w <- data.yau2[9, c(1, 5:177)]
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

# Plot tPCBs --------------------------------------------------------------
# (1) Amanda
selected_values <- c(Air_Amanda = tPCB.conc.air['Conc.Air.Amanda'],
                     WB_Amanda_R = tPCB.conc.wb['wb.amanda.r'],
                     WB_Amanda_L = tPCB.conc.wb['wb.amanda.l'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Amanda R",
                            "Sample 2: WB Amanda L")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.amanda <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (2) Kay
selected_values <- c(Air_Kay = tPCB.conc.air['Conc.Air.Kay'],
                     WB_Kay = tPCB.conc.wb['wb.kay'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Kay")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.kay <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (3) Yau
selected_values <- c(Air_Yau.1 = tPCB.conc.air['Conc.Air.Yau.1st'],
                     WB_Yau_1st = tPCB.conc.wb['wb.yau.1st'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 1st", "Sample 1: WB Yau 1st")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.1 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.2 = tPCB.conc.air['Conc.Air.Yau.2nd'],
                     WB_Yau_1st = tPCB.conc.wb['wb.yau.2nd'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 2nd", "Sample 1: WB Yau 2nd")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.2 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.w = tPCB.conc.air['Conc.Air.Yau.w'],
                     WB_Yau_nw = tPCB.conc.wb['wb.yau.nw'],
                     WB_Yau_w = tPCB.conc.wb['wb.yau.w'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Yau nw",
                            "Sample 2: WB Yau w")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.3 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# Plot all at the same time
tPCB <- grid.arrange(plot.amanda, plot.kay, plot.yau.1, plot.yau.2, plot.yau.3,
             ncol = 3, top = "SumPCBs")

# Save plot in folder
ggsave("Output/Plots/ComparisontPCB.png",
       plot = tPCB, width = 15, height = 10, dpi = 500)

# Select congeners --------------------------------------------------------
tPCB.conc.wb <- rowSums(conc.wb[1:50], na.rm = TRUE)
# Convert list column to numeric, replacing NA values with 0
conc_air_common_numeric <- as.data.frame(sapply(conc_air_common, as.numeric))
tPCB.conc.air <- colSums(conc_air_common_numeric[1:50, ], na.rm = TRUE)

print(tPCB.conc.wb)
print(tPCB.conc.air)

# Plot selected PCBs ------------------------------------------------------
# (1) Amanda
selected_values <- c(Air_Amanda = tPCB.conc.air['Conc.Air.Amanda'],
                     WB_Amanda_R = tPCB.conc.wb['wb.amanda.r'],
                     WB_Amanda_L = tPCB.conc.wb['wb.amanda.l'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Amanda R",
                            "Sample 2: WB Amanda L")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.amanda <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (2) Kay
selected_values <- c(Air_Kay = tPCB.conc.air['Conc.Air.Kay'],
                     WB_Kay = tPCB.conc.wb['wb.kay'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Kay")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.kay <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (3) Yau
selected_values <- c(Air_Yau.1 = tPCB.conc.air['Conc.Air.Yau.1st'],
                     WB_Yau_1st = tPCB.conc.wb['wb.yau.1st'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 1st", "Sample 1: WB Yau 1st")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.1 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.2 = tPCB.conc.air['Conc.Air.Yau.2nd'],
                     WB_Yau_1st = tPCB.conc.wb['wb.yau.2nd'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 2nd", "Sample 1: WB Yau 2nd")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.2 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.w = tPCB.conc.air['Conc.Air.Yau.w'],
                     WB_Yau_nw = tPCB.conc.wb['wb.yau.nw'],
                     WB_Yau_w = tPCB.conc.wb['wb.yau.w'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Yau nw",
                            "Sample 2: WB Yau w")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.3 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# Plot all at the same time
PCBn <- grid.arrange(plot.amanda, plot.kay, plot.yau.1, plot.yau.2, plot.yau.3,
             ncol = 3, top = "Sum 50 PCBs")

# Save plot in folder
ggsave("Output/Plots/ComparisontPCB50.png",
       plot = PCBn, width = 15, height = 10, dpi = 500)

# Concentration calculation using cte SR ----------------------------------
# Divide each element in wb_common cte SR
wb_div_sr <- wb_common/0.35
# Divide further by corresponding value in wb.5d$time.day
conc.wb.cte <- wb_div_sr/wb_time_day
rownames(conc.wb.cte) <- c('wb.amanda.r', 'wb.amanda.l', 'wb.kay', 'wb.yau.1st',
                       'wb.yau.2nd', 'wb.yau.nw', 'wb.yau.w')

# Match congeners in both dataset using cte SR ----------------------------
# Ensure both data frames have matching congener order
common_congener_order <- intersect(names(conc.wb.cte), rownames(conc.air))
# Find indices of matching row names in conc.air
matching_indices <- match(common_congener_order, rownames(conc.air))
# Subset conc.air to include only the rows with matching row names
conc_air_common <- conc.air[matching_indices, ]

# Sum PCBs (2) ------------------------------------------------------------
tPCB.conc.wb.cte <- rowSums(conc.wb.cte, na.rm = TRUE)
# Convert list column to numeric, replacing NA values with 0
conc_air_common_numeric <- as.data.frame(sapply(conc_air_common, as.numeric))
tPCB.conc.air <- colSums(conc_air_common_numeric, na.rm = TRUE)

print(tPCB.conc.wb.cte)
print(tPCB.conc.air)

# Plot tPCBs (2) ----------------------------------------------------------
# (1) Amanda
selected_values <- c(Air_Amanda = tPCB.conc.air['Conc.Air.Amanda'],
                     WB_Amanda_R = tPCB.conc.wb.cte['wb.amanda.r'],
                     WB_Amanda_L = tPCB.conc.wb.cte['wb.amanda.l'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Amanda R cte",
                            "Sample 2: WB Amanda L cte")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.amanda <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (2) Kay
selected_values <- c(Air_Kay = tPCB.conc.air['Conc.Air.Kay'],
                     WB_Kay = tPCB.conc.wb.cte['wb.kay'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Kay cte")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.kay <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# (3) Yau
selected_values <- c(Air_Yau.1 = tPCB.conc.air['Conc.Air.Yau.1st'],
                     WB_Yau_1st = tPCB.conc.wb.cte['wb.yau.1st'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 1st", "Sample 1: WB Yau 1st cte")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.1 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.2 = tPCB.conc.air['Conc.Air.Yau.2nd'],
                     WB_Yau_1st = tPCB.conc.wb.cte['wb.yau.2nd'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration 2nd", "Sample 1: WB Yau 2nd cte")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.2 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

selected_values <- c(Air_Yau.w = tPCB.conc.air['Conc.Air.Yau.w'],
                     WB_Yau_nw = tPCB.conc.wb.cte['wb.yau.nw'],
                     WB_Yau_w = tPCB.conc.wb.cte['wb.yau.w'])

# Rename the categories directly within the vector for clarity
names(selected_values) <- c("Air Concentration", "Sample 1: WB Yau nw cte",
                            "Sample 2: WB Yau w cte")

plot_data <- data.frame(Category = names(selected_values),
                        Value = as.numeric(selected_values))

plot.yau.3 <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab(expression(bold("Air "*Sigma*"PCB (ng/m3)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# Plot all at the same time
tPCB <- grid.arrange(plot.amanda, plot.kay, plot.yau.1, plot.yau.2, plot.yau.3,
                     ncol = 3, top = "Sum25PCBs (SR = 1.5 m3/d)")

# Save plot in folder
ggsave("Output/Plots/Comparison25PCBcte1-5.png",
       plot = tPCB, width = 15, height = 10, dpi = 500)




