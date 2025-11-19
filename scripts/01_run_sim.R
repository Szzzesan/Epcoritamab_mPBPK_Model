library(mrgsolve)
library(tidyverse)

# 1. Load the model from the file
mod <- mread("epc_combined", "model")

# 2. Set up the dose (48 mg SC)
dose_mg <- 48
MW <- 145623 # g/mol
dose_nmol <- (dose_mg / 1000 / MW) * 1e9

# Event 1: Administration
drug_dose <- ev(amt = dose_nmol, cmt = "DEPOT", time = 0)

# Event 2: Injection Effect Trigger
# the paper says this is set to 1 at injection time
inj_signal <- ev(amt = 1, cmt = "INJ", time = 0)

# Combine them into one dosing regimen
total_dose <- drug_dose + inj_signal

# 3. Run the simulation
V_blood_mm3 <- 5e6
total_TC <- 1469 * V_blood_mm3
total_BC <- 236 * V_blood_mm3
# initial state (time 0 - 28 days)
sim_output <- mod %>%
  ev(total_dose) %>%
  # Initialize cell counts at baseline so they don't start at 0
  init(TC_BLOOD = total_TC, BC_BLOOD = total_BC,
       TC_SPLEEN = total_TC * 30, BC_SPLEEN = total_BC * 30) %>% # Rough approx for tissue baseline
  mrgsim(end = 28, delta = 0.1) %>%
  as_tibble()

# 4. Visualization
# # Plot A: Pharmacokinetics (Drug Conc)
# p1 <- ggplot(sim_output, aes(x = time, y = C_plasma)) +
#   geom_line(color = "blue", linewidth = 1) +
#   labs(title = "PK: Plasma Concentration", y = "Conc (nM)", x = "Time (Days)") +
#   theme_bw()
# 
# # Plot B: Lymphocyte Trafficking (T-Cells in Blood)
# # We expect a sharp drop immediately after dosing (The "Injection Effect")
# p2 <- ggplot(sim_output, aes(x = time, y = TC_BLOOD)) +
#   geom_line(color = "red", linewidth = 1) +
#   labs(title = "PD: Blood T-Cell Count",
#        subtitle = "Note the transient drop due to trafficking (INJ effect)",
#        y = "Cells / mm3", x = "Time (Days)") +
#   theme_bw()

# Plot C: Trimer Formation (in blood)
p3 <- ggplot(sim_output, aes(x = time, y = C_Trimer)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "PD: Trimer Concentration",
       subtitle = "The 'Hooks' connecting T-cells and B-cells",
       y = "Conc (nM)", x = "Time (Days)") +
  theme_bw()

# print(p1)
# print(p2)
print(p3)
