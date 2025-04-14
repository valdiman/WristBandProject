## Script to calculate individual PCB "personal" sampling rates
# for silicone wristbands.
# nd non-dominant hand
# d dominant hand

# Install packages
install.packages("readxl")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("deSolve")
install.packages("gridExtra")
install.packages("minpack.cl")

# Load libraries
{
  library(readxl)
  library(reshape2) # organize data
  library(ggplot2)
  library(deSolve) # solving differential equations
  library(gridExtra)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
  library(rlang)
  library(minpack.lm)
}

# Read data from excel ----------------------------------------------------
data.amanda <- data.frame(read_excel("Data/Amanda.xlsx", sheet = "Sheet1",
                             col_names = TRUE, col_types = NULL))
data.kay <- data.frame(read_excel("Data/Kay.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet1",
                                     col_names = TRUE, col_types = NULL))
data.yau2 <- data.frame(read_excel("Data/Yau.xlsx", sheet = "Sheet2",
                                  col_names = TRUE, col_types = NULL))
logKoa <- data.frame(read_excel("Data/logKoa.xlsx", sheet = "logKoa",
                                col_names = TRUE, col_types = NULL))

# Calculate air concentration and read Amanda wb data ---------------------------------
# WBs were used to calculate PCB concentration
# Both hands (d and nd)
# triplicates for 4.3 days were deployed
# sampling rate of 0.5 m3/d was used for static WBs
{
  # Select WBs to calculate air concentration
  data.amanda.1 <- data.amanda[1:3, ]
  # Average 3 WBs. NA values not included in the calculations
  data.amanda.2 <- colMeans(data.amanda.1[, 3:175], na.rm = TRUE)
  # Calculate air concentration in ng/m3
  # = massWB/(0.5*time.day)
  conc <- data.amanda.2/(0.5*data.amanda[1,1])
  # Amanda right wb data
  a.wb.r <- data.amanda[4:8, 3:175]
  # Amanda left wb data
  a.wb.l <- data.amanda[9:13, 3:175]
}

time <- data.amanda$time.day[4:8] * 24 # [h]
pcb.ind <- "PCB8"

cair.pcbi <- conc[[pcb.ind]] # [ng/m3]
a.wb.r.pcbi <- a.wb.r[[pcb.ind]] # [ng/WB g]

# WB physical parameters
Vwb <- 4.73 * 10^-6 # [m3]
mwb <- 4 # [g]
dwb <- 1080000 # [g/m3]

# Define ODE with constant Cair
predict_Xwb <- function(pars, times, cair.pcbi, dwb) {
  ku <- pars["ku"]
  ke <- pars["ke"]
  
  ode_func <- function(t, state, parameters) {
    Xwb <- state[1]
    dXwb <- ku / dwb * cair.pcbi - ke * Xwb
    return(list(c(dXwb)))
  }
  
  state0 <- c(Xwb = 0)  # Xenm @ time 0 = 0
  
  ode_out <- ode(y = state0, times = times, func = ode_func, parms = NULL)
  Xwb_pred <- ode_out[, "Xwb"]
  return(Xwb_pred)
}

# Residuals function for fitting
fit_model <- function(pars, times, Xwb_obs, cair.pcbi, dwb) {
  sim_Xwb <- predict_Xwb(pars, times, cair.pcbi, dwb)
  return(sim_Xwb - Xwb_obs)
}

# Starting parameters
start_pars <- c(ku = 10000, ke = 0.1) # [1/h]

# Fit using constant Cair
fit <- nls.lm(
  par = start_pars,
  fn = fit_model,
  lower = c(0, 0),
  control = nls.lm.control(maxiter = 200),
  times = time,
  Xwb_obs = a.wb.r.pcbi,
  cair.pcbi = cair.pcbi,
  dwb = dwb
)

summary(fit)

# Predicted values at original time points
Xwb_fitted <- predict_Xwb(fit$par, time, cair.pcbi, dwb)

# Model evaluation metrics (performance metrics)
# Root mean squared error (RMSE)
RMSE <- sqrt(mean((a.wb.r.pcbi - Xwb_fitted)^2))
# Residual sum of squares (RSS)
RSS <- sum((a.wb.r.pcbi - Xwb_fitted)^2) # [ng/g]
# Total sum of squares (TSS)
TSS <- sum((a.wb.r.pcbi - mean(a.wb.r.pcbi))^2)
# R-squared
R2 <- 1 - RSS / TSS

# Sampler parameters
# Kenm calculations
ku <- fit$par[1]
ke <- fit$par[2]
logKwb <- log10(ku / ke * 1000^2 / dwb) # [m3/m3]
# Sampling rate
Rs <- ku * Vwb * 24 # [m3/d]
# Time to reach 90% equilibrium
t90 <- log(10) / ke #[h]

# Create a summary data frame
model_summary <- data.frame(
  pcb.ind = pcb.ind,
  ku = ku,
  ke = ke,
  logKenm = logKwb,
  Rs = Rs,
  t90 = t90,
  RMSE = RMSE,
  R2 = R2,
  row.names = NULL
)

# Add units to column names
colnames(model_summary) <- c(
  "congener",
  "ku (1/h)",
  "ke (1/h)",
  "logKwb (L/kg)",
  "Rs (m3/d)",
  "t90 (h)",
  "RMSE",
  "R2"
)

# Print the summary for verification
print(model_summary)

# Predict for plotting
pars_fit <- fit$par
time_smooth <- seq(min(time), max(time), by = 0.1)
Xwb_pred <- predict_Xwb(pars_fit, time_smooth, cair.pcbi, dwb)

# Plot
obs_df <- data.frame(time = time, Observed = a.wb.r.pcbi)
smooth_df <- data.frame(time = time_smooth, Predicted = Xwb_pred)

plot.uptake <- ggplot() +
  geom_point(data = obs_df, aes(x = time, y = Observed), shape  = 21,
             color = "black", size = 2.5) +
  geom_line(data = smooth_df, aes(x = time, y = Predicted),
            color = "black", linewidth = 0.4) +
  theme_bw() +
  labs(x = expression(bold("Time (hours)")),
       y = bquote(bold("ng" ~.(pcb.ind) ~ "/g WB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10, vjust = 0.1),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10),
        aspect.ratio = 1,
        panel.grid = element_blank())

plot.uptake



