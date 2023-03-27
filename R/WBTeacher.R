# R code to analyze school teacher wristband data

# Install packages
install.packages("reshape2")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("readxl") #say no!

# Library
{
  library(readxl)
  library(ggplot2)
  library(reshape2) # For melt function
  library(tidyverse) # Data manipulation 
  library(dplyr) # performs %>%
}

# Read data ---------------------------------------------------------------
# Read blank data from excel
bl <- data.frame(read_excel("Data/BlanksWBT.xlsx", sheet = "blanks",
                             col_names = TRUE, col_types = NULL))
# Remove entire "bad/yellow" transition
s.1 <- data.frame(read_excel("Data/SamplesWBT.xlsx", sheet = "samples",
                            col_names = TRUE, col_types = NULL))
s.2 <- data.frame(read_excel("Data/SamplesWBT.xlsx", sheet = "samplesv02",
                             col_names = TRUE, col_types = NULL))

# Distribution analysis ---------------------------------------------------
# Remove metadata from blank data
bl.1 <- subset(bl, select = -c(Congener.Blank))
# Look at the distribution of the blank data
# Create matrix to storage data
normality <- matrix(nrow = length(bl.1[1,]), ncol = 2)

for (i in 1:length(bl.1[1, ])) {
  normality[i, 1] <- shapiro.test(bl.1[,i])$p.value
  normality[i, 2] <- shapiro.test(log10(bl.1[,i]))$p.value
}
  
# Just 3 significant figures
normality <- formatC(signif(normality, digits = 3))
# Add congener names
congeners <- colnames(bl.1)
# Include congener names
normality <- cbind(congeners, normality)
# Change column names
colnames(normality) <- c("Congener", "shapiro.normal", "shapiro.log10")
# Create Q-Q plot for individual PCB congeners
{
  qqnorm(bl.1$PCB9, main = "Concentration (ng/g)")
  qqline(bl.1$PCB9)
}

# Bar plot to visualize Shapiro test (normality)
# Organize data to be plotted
norma <- data.frame(normality)
norma$Congener <- as.character(norma$Congener)
norma$shapiro.normal<- as.numeric(as.character(norma$shapiro.normal))
norma$shapiro.log10 <- as.numeric(as.character(norma$shapiro.log10))

# Plot
# Values less than -log10(0.05) are not significant
# Normal data
ggplot(norma, aes(x = Congener, y = -log10(shapiro.normal))) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "black") +
  theme_test() +
  theme(aspect.ratio = 2/12) +
  ylab("-log10 (p-value)") +
  theme(axis.title.y = element_text(face = "bold", size = 8)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "red",
             linewidth = 0.8)

# log10 data
ggplot(norma, aes(x = Congener, y = -log10(shapiro.log10))) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "black") +
  theme_test() +
  theme(aspect.ratio = 2/12) +
  ylab("-log10 (p-value)") +
  theme(axis.title.y = element_text(face = "bold", size = 8)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = -log10(0.05), color = "red",
             linewidth = 0.8)
  
# Export data
write.csv(normality, file = "Output/Data/csv/normality.csv")

# Calculate LOQ -----------------------------------------------------------
# Create LOQ, i.e., upper 95 CI% (mean + 1.96*sd/(n)^0.5)
loq <- colMeans(bl.1) + 1.96*sapply(bl.1, sd)/sqrt(20)
loq <- data.frame(t(loq))

# Make comparison between s.1 and loq
# If s.1 > loq, then s.1, if not 0
# Create matrix to storage s.1 or loq values in s.2
 <- matrix(NA, nrow = dim(s.1)[1], ncol = dim(s.1)[2])
# Do comparison
for(i in 1:dim(s.1)[1]) {
  for(j in 1:dim(s.1)[2]) {
    if (log10(s.1[i, j]) > loq[j]) {
      s.2[i, j] <- s.1[i, j]
    } else {
      s.2[i, j] <- 0
    }
  }
}



