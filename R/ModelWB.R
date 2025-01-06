# Model WB

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")
install.packages("gridExtra")

# Load libraries
{
  library(dplyr) # organize data
  library(reshape2) # organize data
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
  library(gridExtra)
}

# Uptake WB function ------------------------------------------------------
uptakeWB.PCB52 = function(t, state, parms){
  
  Kwb <- 10^4
  Vwb <- 1
  
  # Sampling rate
  sr <- parms$sr
  
  # concentrations
  Mwb <- state[1]
  Vef <- state[2]
  
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  dVefdt <- sr / Vwb * (Vwb - Vef / Kwb)
  
  # The computed derivatives are returned as a list
  return(list(c(dMwbdt, dVefdt)))
}

# Initial conditions and run function
Ca <- 10 # [ng/m3]

cinit <- c(Mwb = 0, Vef = 0) # [ng/m3]
parms <- list(sr = 1)
t <- seq(0, 7, 0.1) # [d]
# Run the ODE function without specifying parms
model <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms)
head(model)

# Transform model to data.frame and get air concentation from Mwb
model.2 <- as.data.frame(model)
colnames(model.2) <- c("time", "Mwb", "Vef")
model.2$Cwb <- model.2$Mwb / model.2$Vef



