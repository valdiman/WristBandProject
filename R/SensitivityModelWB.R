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
Ca_values <- seq(0.1, 10, by = 2)  # Vary Ca from 0.1 to 10
sr_values <- seq(0.1, 2, by = 0.4) # Vary sr from 0.1 to 2

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
  # Check if result has required columns
  if ("time" %in% colnames(res) && "Mwb" %in% colnames(res)) {
    data.frame(
      time = res[, "time"],
      Mwb = res[, "Mwb"],  # Keep Mwb values
      scenario = scenario
    )
  } else {
    NULL  # Skip invalid results
  }
}))

# Plot the sensitivity results
ggplot(sensitivity_df, aes(x = time, y = Mwb, color = scenario)) +
  geom_line() +
  theme_bw() +
  labs(title = "Sensitivity Analysis of Mwb",
       x = "Time (days)", y = "Mwb (ng)")

# Extract the values of sr and Ca from the scenario name
sensitivity_df$sr <- as.numeric(sub(".*sr([0-9.]+)_.*", "\\1", sensitivity_df$scenario))
sensitivity_df$Ca <- as.numeric(sub("Ca([0-9.]+)_.*", "\\1", sensitivity_df$scenario))

# Correlation between Mwb and sr, and between Mwb and Ca
cor_sr <- cor(sensitivity_df$Mwb, sensitivity_df$sr)
cor_Ca <- cor(sensitivity_df$Mwb, sensitivity_df$Ca)

# Print correlation results
print(paste("Correlation between Mwb and sr:", cor_sr))
print(paste("Correlation between Mwb and Ca:", cor_Ca))

