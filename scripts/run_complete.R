library(mrgsolve)
library(tidyverse)

# 1. Load the Master Model
# Ensure the file is named "epc_complete.cpp" inside the "model" folder
mod <- mread("epc_complete", "model", preclean = TRUE)

# 2. Set up Total Cell Numbers & Define Dosing (48 mg SC)
# 48 mg -> nmol conversion (mg -> g -> mol -> nmol)
dose_nmol <- (48 / 1000 / 145623) * 1e9 

# Combined Dose: Drug (DEPOT) + Injection Effect (INJ)
total_dose <- ev(amt = dose_nmol, cmt = "DEPOT", time = 0) + 
  ev(amt = 1, cmt = "INJ", time = 0)

V_blood_mm3 <- 5e6 
total_TC <- 1469 * V_blood_mm3
total_BC <- 236  * V_blood_mm3

# 3. Run Simulation
# Note: We use a high maxsteps because binding kinetics are fast
sim <- mod %>%
  ev(total_dose) %>%
  param(scale_binding = 1) %>%
  init(
    TC_BLOOD = total_TC, 
    BC_BLOOD = total_BC, 
    TC_SPLEEN = total_TC * 30,
    BC_SPLEEN = total_BC * 30, 
    TC_NODE = total_TC * 0.03, 
    BC_NODE = total_BC * 0.03, 
    TC_LYMPH = total_TC * 19, 
    BC_LYMPH = total_BC * 19, 
    AF_TC = 1, 
    AF_BC = 1, 
    INJ = 0
  ) %>%
  # Run for 12 weeks (84 days)
  mrgsim(end = 84, delta = 0.1, 
         atol = 1e-3, rtol = 1e-3,
         maxsteps = 1e9) %>%
  as_tibble()
print('Finish running the simulation')

# --- 4. DEBUGGING PLOTS (The Chain of Causality) ---

# A. PK Check: Is drug getting into the system?
p1 <- ggplot(sim, aes(x = time, y = C_plasma)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "1. PK Check: Plasma Conc", y = "Conc (nM)", x = "Time (Days)") +
  theme_bw()
print(p1)

# B. Trafficking Check: Are T-Cells leaving the blood?
p2 <- ggplot(sim, aes(x = time, y = TC_BLOOD)) +
  geom_line(color = "orange", linewidth = 1) +
  labs(
    title = "2. Trafficking Check: Blood T-Cells", 
    subtitle = "Note the immediate drop (Injection Effect) and slow recovery (Homeostasis)",
    y = "Cell Count", 
    x = "Time (Days)"
  ) +
  theme_bw()
print(p2)

# C. Binding Check: Are Trimers forming in the blood?
p3 <- ggplot(sim, aes(x = time, y = C_Trimer_BL_nM)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "3. Binding Check: Blood Trimer", y = "Conc (nM)", x = "Time (Days)") +
  theme_bw()
print(p3)

p3_zoomed <- ggplot(sim, aes(x = time, y = C_Trimer_BL_nM)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "3. Binding Check: Blood Trimer (Zoomed in)",
       subtitle = "The second smaller hump is likely driven by T-cell expansion",
       y = "Conc (nM)", x = "Time (Days)") +
  coord_cartesian(ylim = c(0, 1e-6)) +
  theme_bw()
print(p3_zoomed)

p_bcell <- ggplot(sim, aes(x = time, y = BC_BLOOD)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    title = "Diagnostic 1: B-Cell Depletion",
    subtitle = "B-Cells get depleted as an intended therapeutic goal of Epc",
    y = "B-Cell Count (Blood)",
    x = "Time (Days)"
  ) +
  scale_y_continuous(labels = scales::scientific) +
  # coord_cartesian(xlim = c(0, 10)) + # Zoom in on the crash
  theme_bw()
print(p_bcell)


# D. Activation Check: Are T-cells expanding?
p4 <- ggplot(sim, aes(x = time, y = Total_pATC)) +
  geom_line(color = "purple", linewidth = 1) +
  labs(title = "4. Activation Check: Total Activated T-Cells", y = "Cell Count", x = "Time (Days)") +
  theme_bw()
print(p4)

# E. Efficacy Check: Is the Tumor shrinking?
p5 <- ggplot(sim, aes(x = time, y = Tumor_Vol_cm3)) +
  geom_line(color = "red", linewidth = 1) +
  labs(title = "5. Efficacy: Tumor Volume", y = "Volume (cm3)", x = "Time (Days)") +
  theme_bw()
print(p5)

