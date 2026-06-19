# Script used to model wristband uptake of PCB52 under variable airborne concentrations.
# The results are not included in the manuscript.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("ggplot2")
install.packages("deSolve")
install.packages("tidyr")

# Load libraries
{
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
  library(tidyr)
}

# (i) Constant sr and constant Ca -----------------------------------------
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
  geom_line() +
  geom_point(data = obs, aes(x = time.obs, y = vef.obs), size  = 4) +
  theme_bw() +
  labs(x = "time (days)", y = "Effective volume PCB 52 (m3)") +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14))

# Same as (i) but WB rotating
uptakeWBr.PCB52 = function(t, state, parms){
  
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
Ca <- 33 # [ng/m3]

cinit <- c(Mwb = 0, Vef = 0) # [ng/m3]
parms <- list(sr = 1.1) # from SamplingRates.R
t <- seq(0, 5, 0.1) # [d]
# Run the ODE function without specifying parms
model.r <- ode(y = cinit, times = t, func = uptakeWBr.PCB52, parms = parms)
head(model.r)

# Transform model to data.frame and get air concentation from Mwb
model.r.df <- as.data.frame(model.r)
colnames(model.r.df) <- c("time", "Mwb", "Vef")
model.r.df$Cwb <- model.r.df$Mwb / model.r.df$Vef

# Include observations from uptake studies
vef.r.obs <- c(1.188, 2.291, 3.508, 4.08) # [m3] Vef
time.r.obs <- c(1.025, 2.03, 3.053, 4.0625) # [d] time
obs.r <- cbind(vef.r.obs, time.r.obs)

# Plot Vef vs. time
ggplot(data = model.r.df, aes(x = time, y = Vef)) +
  geom_line() +
  geom_point(data = obs.r, aes(x = time.r.obs, y = vef.r.obs), size = 4) +
  theme_bw() +
  labs(x = "time (days)", y = "Effective volume PCB 52 (m3)") +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14))

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
  
  Kwb <- 10^7.4 # [m3wb/m3air] from 4.
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Variables
  Mwb <- state[1] # [ng]
  Vef <- state[2] # [m3]
  
  # Calculate the current hour (time in hours, 0 to 23)
  hour <- t %% 24 # Get the current hour (0 to 23)
  
  # Apply the conditions for sr
  sr <- ifelse(hour < 7, 0.5,                  # 0:00 to 7:00
               ifelse(hour < 22, 0.9,         # 7:00 to 22:00
                      0.5))                   # 22:00 to 24:00
  
  # Apply the conditions for Ca
  Ca <- ifelse(hour < 8, 0.2,                            # 0:00 to 8:00 home
               ifelse(hour < 8.5, 0.0227,                # 8:00 to 8:30 outdoor/commute
                      ifelse(hour < 17, 1.4,             # 8:30 to 17:00 office
                             ifelse(hour < 17.5, 0.0227, # 17:00 to 17:30 outdoor/commute
                                    0.2))))              # 17:30 to 24:00 home
  
  
  # Differential equations
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  dVefdt <- sr / Vwb * (Vwb - Vef / Kwb)  # [m3/d]
  
  # Return computed derivatives as a list, including Ca and sr
  return(list(c(dMwbdt, dVefdt), Ca = Ca, sr = sr))
}

# Initial conditions
cinit <- c(Mwb = 0, Vef = 0) # Initial values for Mwb and Vef

# Time sequence for 5 days (starting at 8 AM on Monday and ending at 5 PM on Friday)
start_time <- 9   # Starting at 9 AM on Monday (first day)
end_time <- 4 * 24 + 17  # Ending at 5 PM on Friday (5th day)

# Generate time sequence with 0.25 hour steps
t_seq <- seq(start_time, end_time, by = 0.25)

# Solve the differential equations using ode
model.3 <- ode(y = cinit, times = t_seq, func = uptakeWB3.PCB52, parms = NULL)

# View the result
head(model.3)
# Transform model to data.frame
model.df.3 <- as.data.frame(model.3)
# Add Cwb
model.df.3$Cwb <- model.df.3$Mwb / model.df.3$Vef
# To check hours of the day
model.df.3$hour_of_day <- model.df.3$time %% 24

# Plot Vef vs. time
ggplot(data = model.df.3, aes(x = time, y = Vef)) +
  geom_line(linetype = "dashed") +
  theme_classic()

# Calculate the average Ca over the 5 days
Ca_avg <- mean(model.df.3$Ca)

# (1) Cwb, Ca vs. time (observations not included)
# Start time: Monday, January 6th, 2025, at 09:00 AM
start_time <- as.POSIXct("2025-01-06 09:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Generate a sequence of times with an increment of 0.25 hours (15 minutes)
model.df.3$new_datetime <- start_time + (0:416) * 3600 * 0.25  # 0.25 hours = 15 minutes

# Format the new datetime to display only the day and hour (without minutes/seconds)
model.df.3$new_formatted_datetime <- format(model.df.3$new_datetime, "%a %H:%M")

# Subset the data to include only every 24th row
model.df.3_subset <- model.df.3[seq(1, nrow(model.df.3), by = 32), ]

# Plot the data
plotPCB52Model <- ggplot(data = model.df.3, aes(x = time)) +
  geom_line(aes(y = Cwb, color = "Cwb"), linetype = "dashed", linewidth  = 1) +
  geom_line(aes(y = Ca, color = "Ca")) +
  geom_line(aes(y = Ca_avg, color = "Ca_avg"), linetype = "dotted", linewidth = 1) +
  scale_color_manual(name = "",
                     values = c("Cwb" = "blue", "Ca" = "#009E73", "Ca_avg" = "black"),
                     labels = c("Cwb" = "Estimate air concentration from full-day WB", 
                                "Ca" = "Air concentration",
                                "Ca_avg" = "Average air concentration")) +
  scale_x_continuous(name = "Time (hours)", 
                     breaks = seq(0, 142, by = 12),  # Customize breaks for the primary axis
                     labels = seq(0, 142, by = 12)) + 
  geom_text(data = model.df.3_subset, aes(x = time, y = 1.6, label = new_formatted_datetime), 
            angle = 90, hjust = 1, vjust = 0.5, size = 2.5) +
  labs(y = "PCB 52 (ng/m3)") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
  legend.text = element_text(size = 8, face = "bold"))

# Print the plot
print(plotPCB52Model)

# Save plot in folder
ggsave("Output/Plots/Model/PCB52Model.png",
       plot = plotPCB52Model, width = 8, height = 5, dpi = 500)

# Add Cwb from volunteers
Cwb.obs <- c(0.58, 0.6, 2.62, 1.27, 0.88)
time.obs <- c(3.2 * 24, 3.5 * 24, 3.3 * 24, 3.8* 24, 3.6 * 24) # [d] time
obs <- cbind(Cwb.obs, time.obs)
# Mi, Ea, Ya, An & Xu
rownames(obs) <- c("Vol. 1", "Vol. 2", "Vol. 3", "Vol. 4", "Vol. 5")

# (2) Cwb, Ca vs. time
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
  scale_x_continuous(name = "Time (hours)", 
                     breaks = seq(0, 142, by = 12), # Customize break points here
                     labels = seq(0, 142, by = 12)) +  # Customize labels (optional)
  labs(y = "PCB 52 (ng/m3)") +
  theme(axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14))

