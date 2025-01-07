# Model WB for PCB 52

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

# (i) Constant sr and constant Ca -----------------------------------------
# Uptake WB function ------------------------------------------------------
uptakeWB.PCB52 = function(t, state, parms){
  
  Kwb <- 10^7.4 # [m3wb/m3air] from Frederiksen 2022
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Sampling rate
  sr <- parms$sr # [m3/d]
  
  # Variables
  Mwb <- state[1] # [ng]
  Vef <- state[2] # [m3]
  
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  dVefdt <- sr / Vwb * (Vwb - Vef / Kwb) # [m3/d]
  
  # The computed derivatives are returned as a list
  return(list(c(dMwbdt, dVefdt)))
}

# Initial conditions and run function
Ca <- 20 # [ng/m3]

cinit <- c(Mwb = 0, Vef = 0) # [ng/m3]
parms <- list(sr = 0.53)
t <- seq(0, 5, 0.1) # [d]
# Run the ODE function without specifying parms
model <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms)
head(model)

# Transform model to data.frame and get air concentation from Mwb
model.df <- as.data.frame(model)
colnames(model.df) <- c("time", "Mwb", "Vef")
model.df$Cwb <- model.df$Mwb / model.df$Vef

# Include observations from uptake studies
vef.obs <- c(0.621, 1.114, 1.66, 1.924, 2.7) # [m3] Vef
time.obs <- c(1, 2, 2.97, 3.968, 4.982) # [d] time
obs <- cbind(vef.obs, time.obs)

# Plot Vef vs. time
ggplot(data = model.df, aes(x = time, y = Vef)) +
  geom_line(linetype = "dashed") +
  geom_point(data = obs, aes(x = time.obs, y = vef.obs)) +
  theme_classic()
  
# (ii) Variable sr and constant Ca -----------------------------------------
uptakeWB2.PCB52 <- function(t, state, parms) {
  
  Kwb <- 10^7.4 # [m3wb/m3air] from Frederiksen 2022
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Time-varying sampling rate (sr)
  sr <- ifelse((t %% 1) < (16 / 24), 0.9, 0.5) # 1 for first 16 hours, 0.5 for next 8 hours
  
  # Variables
  Mwb <- state[1] # [ng]
  Vef <- state[2] # [m3]
  
  # Differential equations
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  dVefdt <- sr / Vwb * (Vwb - Vef / Kwb)  # [m3/d]
  
  # Return computed derivatives as a list
  return(list(c(dMwbdt, dVefdt)))
}

# Initial conditions
Ca <- 20 # [ng/m3]
cinit <- c(Mwb = 0, Vef = 0) # Initial values for Mwb and Vef
# Time sequence for 5 days (in hours)
t <- seq(0, 5, 0.1) # 5 days in hours with a 0.1-hour step
# Run the ODE solver
model.2 <- ode(y = cinit, times = t, func = uptakeWB2.PCB52, parms = NULL)
head(model.2)

# Transform model to data.frame and get air concentation from Mwb
model.df.2 <- as.data.frame(model.2)
colnames(model.df.2) <- c("time", "Mwb", "Vef")
model.df.2$Cwb <- model.df.2$Mwb / model.df.2$Vef

# Plot Vef vs. time
ggplot(data = model.df.2, aes(x = time, y = Vef)) +
  geom_line(linetype = "dashed") +
  theme_classic()

# (iii) Variable sr and Ca ------------------------------------------------
# Ca from BSB student room
uptakeWB3.PCB52 <- function(t, state, parms) {
  
  Kwb <- 10^7.4 # [m3wb/m3air] from Frederiksen 2022
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Time-varying sampling rate (sr)
  sr <- ifelse((t %% 1) < (16 / 24), 0.9, 0.5) # 0.9 for first 16 hours, 0.5 for next 8 hours
  
  # Time-varying air concentration (Ca)
  Ca <- ifelse((t %% 1) < (8 / 24), 1.4, # 1.4 for 8 hours at work
               ifelse((t %% 1) < (9 / 24), 0.0227, # 0.0227 for 1 hour commuting (outdoor)
                      0.1)) # 0.1 for remaining time at home
  
  # Variables
  Mwb <- state[1] # [ng]
  Vef <- state[2] # [m3]
  
  # Differential equations
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  dVefdt <- sr / Vwb * (Vwb - Vef / Kwb)  # [m3/d]
  
  # Return computed derivatives as a list
  return(list(c(dMwbdt, dVefdt), Ca = Ca, sr = sr))
}

# Initial conditions
cinit <- c(Mwb = 0, Vef = 0) # Initial values for Mwb and Vef
# Time sequence for 5 days (in days)
t <- seq(0, 5, 0.01) # 5 days in days with 0.1-day steps
# Run the ODE solver
model.3 <- ode(y = cinit, times = t, func = uptakeWB3.PCB52, parms = NULL)
head(model.3)

# Transform model to data.frame
model.df.3 <- as.data.frame(model.3)
colnames(model.df.3) <- c("time", "Mwb", "Vef", "Ca", "sr")
model.df.3$Cwb <- model.df.3$Mwb / model.df.3$Vef

# Plot Vef vs. time
ggplot(data = model.df.3, aes(x = time, y = Vef)) +
  geom_line(linetype = "dashed") +
  theme_classic()

# Calculate the average Ca over the 5 days
Ca_avg <- mean(model.df.3$Ca)

# Add Cwb from volunteers
Cwb.obs <- c(0.58, 0.6, 2.62, 1.27, 0.88)
time.obs <- c(3.2, 3.5, 3.3, 3.8, 3.6) # [d] time
obs <- cbind(Cwb_obs, time.obs)
# Mi, Ea, Ya, An & Xu
rownames(obs) <- c("Vol. 1", "Vol. 2", "Vol. 3", "Vol. 4", "Vol. 5")

# Cwb, Ca vs. time
ggplot(data = model.df.3, aes(x = time)) +
  geom_line(aes(y = Cwb, color = "Cwb")) +
  geom_line(aes(y = Ca, color = "Ca")) +
  geom_line(aes(y = Ca_avg, color = "Ca_avg"), linetype = "dashed", linewidth = 1) + # Average Ca line
  scale_color_manual(name = "Variables",
                     values = c("Cwb" = "blue", "Ca" = "#009E73", "Ca_avg" = "black"),
                     labels = c("Cwb" = "Estimate Air concentration from WB", 
                                "Ca" = "Air concentration",
                                "Ca_avg" = "Average Air concentration")) +
  theme_bw() +
  geom_point(data = obs, aes(x = time.obs, y = Cwb.obs)) +
  geom_text(data = obs, aes(x = time.obs, y = Cwb.obs, label = rownames(obs)), 
            vjust = -1, color = "black", size = 4) +  # Adding row names as labels
  labs(x = "Time (days)", y = "PCB 52 (ng/m3)") +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14))





