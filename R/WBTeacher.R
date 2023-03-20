# R code to analyze school teacher wristband data

# Install packages
install.packages("userfriendlyscience")
install.packages("reshape2")
install.packages("tidyverse")
install.packages("readxl") #say no!

# Library
{
  library(readxl)
  library(userfriendlyscience) # Perform Tukey test
  library(reshape2) # For melt function
  library(tidyverse) # Data manipulation 
  library(dplyr) # performs %>%
}

# Read data ---------------------------------------------------------------
# Read data from excel
bl <- data.frame(read_excel("Data/BlanksWBT.xlsx", sheet = "blanks",
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
# Export data
write.csv(normality, file = "Output/Data/csv/normality.csv")

# Calculate LOQ -----------------------------------------------------------
# Create LOQ, i.e., upper 95 CI% (mean + 1.96*sd/(n)^0.5)
loq <- colMeans(bl.1) + 1.96*sapply(bl.1, sd)/sqrt(20)
loq <- data.frame(t(loq))


# Make comparison between s.1 and loq
# If s.1 > loq, then s.1, if not 0
# Create matrix to storage s.1 or loq values in s.2
s.2 <- matrix(NA, nrow = dim(s.1)[1], ncol = dim(s.1)[2])
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



