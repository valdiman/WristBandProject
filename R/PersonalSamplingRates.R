## Code to individual PCB calculate "personal" sampling rates
# for silicone wristbands. Sampling rates were calculated with

# Install packages
install.packages("readxl") #say no!
install.packages("gridExtra")
install.packages("ggplot2")

# Load libraries
{
  library(readxl)
  library(ggplot2)
  library(gridExtra)
}

# Read data from excel ----------------------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                             col_names = TRUE, col_types = NULL))
data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))

# Calculate personal sampling rate Amanda ---------------------------------
# Select WBs to calculate air concentration
data.amanda.1 <- data.amanda[1:3,]



# Remove medata
data.amanda.2 <- subset(data.amanda.1, select = -c(time.day:congeners))


