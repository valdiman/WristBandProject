
# Packages and libraries --------------------------------------------------
# Install packages
install.packages("ggplot2")
install.packages("deSolve")

# Load libraries
{
  library(ggplot2) # plotting
  library(deSolve) # solving differential equations
}

# Uptake WB function ------------------------------------------------------
uptakeWB.PCB52 <- function(t, state, parms) {
  
  Kwb <- 10^7.4 # [m3wb/m3air] from Frederiksen 2022
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Extract parameters
  sr <- parms$sr # [m3/d]
  Ca <- parms$Ca # [ng/m3]
  
  # Variables
  Mwb <- state[1] # [ng]
  
  # Define the differential equation for Mwb
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  
  # Return the computed derivative for Mwb
  return(list(c(dMwbdt)))
}

# Define ranges for Ca and sr
Ca_values <- seq(0.1, 25, by = 4.5)  # Vary Ca from 0.1 to 10
sr_values <- seq(0.1, 2, by = 0.35) # Vary sr from 0.1 to 2

# Initialize storage for results
results <- list()

# Initial conditions and time sequence
cinit <- c(Mwb = 0) # [ng/m3], only one state variable (Mwb)
t <- seq(0, 5, 0.1) # [d]

# Loop over Ca and sr values
for (Ca in Ca_values) {
  for (sr in sr_values) {
    parms <- list(sr = sr, Ca = Ca)
    model <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms)
    results[[paste("Ca", Ca, "sr", sr, sep = "_")]] <- model
  }
}

# Extract and summarize results
sensitivity_df <- do.call(rbind, lapply(names(results), function(scenario) {
  res <- results[[scenario]]
  
  # Extract Ca and sr directly from the scenario string
  Ca_value <- as.numeric(gsub("Ca_([0-9.]+)_sr_.*", "\\1", scenario))
  sr_value <- as.numeric(gsub(".*_sr_([0-9.]+)", "\\1", scenario))
  
  # Check if result has required columns
  if ("time" %in% colnames(res) && "Mwb" %in% colnames(res)) {
    data.frame(
      time = res[, "time"],
      Mwb = res[, "Mwb"],  # Keep Mwb values
      Ca = Ca_value,
      sr = sr_value       # Add Ca and sr explicitly
    )
  } else {
    NULL  # Skip invalid results
  }
}))

# Plot Mwb over time, faceted by Ca and colored by sr
ggplot(sensitivity_df, aes(x = time, y = Mwb, color = factor(sr))) +
  geom_line() +
  facet_wrap(~ Ca, scales = "free_y") +
  labs(
    title = "Mwb ~ f(Ca, sr)",
    x = "Time (days)",
    y = "Mwb (ng)",
    color = "sr (m3/d)"
  ) +
  theme_minimal()

# Change Ca and sr at one time.

uptakeWB.PCB52 = function(t, state, parms){
  
  Kwb <- 10^7.4 # [m3wb/m3air] from Frederiksen 2022
  Vwb <- 4.73 * 10^-6 # [m3]. Debossed adult size (www.24hourwristbands.com)
  
  # Sampling rate
  sr <- parms$sr # [m3/d]
  
  # Variables
  Mwb <- state[1] # [ng]
  
  dMwbdt <- sr * (Ca - Mwb / (Kwb * Vwb)) # [ng/d]
  
  # The computed derivatives are returned as a list
  return(list(c(dMwbdt)))
}

# Initial conditions and run function
Ca <- 20 # [ng/m3]
cinit <- c(Mwb = 0) # [ng]
parms <- list(sr = 0.5)
t <- seq(0, 5, 0.01) # [d]
# Baseline
model.bl <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms) # bl = baseline
head(model.bl)
model.bl <- as.data.frame(model.bl)
# sr 50% +
parms <- list(sr = 0.5 * 1.5)
model.sr <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms)
head(model.sr)
model.sr <- as.data.frame(model.sr)
percent_change_sr <- (model.sr$Mwb - model.bl$Mwb) / model.bl$Mwb * 100
percent_change_sr
# Ca 50% +
Ca <- 20 * 1.5 # [ng/m3]
parms <- list(sr = 0.5)
model.Ca <- ode(y = cinit, times = t, func = uptakeWB.PCB52, parms = parms) # bl = baseline
head(model.Ca)
model.Ca <- as.data.frame(model.Ca)
percent_change_Ca <- (model.Ca$Mwb - model.bl$Mwb) / model.bl$Mwb * 100
percent_change_Ca

# Outcome: Both increase Mbw at a very similar amount, 50% more.
