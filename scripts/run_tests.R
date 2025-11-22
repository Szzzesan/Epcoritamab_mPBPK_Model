library(mrgsolve)
library(tidyverse)

# --- 1. SETUP ---
# Calculate Total Cell Counts for Initialization
V_blood_mm3 <- 5e6 
total_TC <- 1469 * V_blood_mm3
total_BC <- 236  * V_blood_mm3

# Ratios from Paper (approximate):
# Blood = 1x (Reference)
# Spleen = 30x (60% vs 2%)
# Node + Lymph = ~19x (38% vs 2%) -> Biologically consistent values would be 0.03 and 19
# At steady state, the flow-in must equal the flow-out. From node to lymph, 
# rate k_{nl} = 1672, so to match the flow: 1672 x Node = 50 x Blood, so node 
# should be 50/1672 = 0.03 x Blood

# Dose: 48 mg
dose_nmol <- (48 / 1000 / 145623) * 1e9

# separate the actual drug dose and 
# the body's reaction to the injection (the injection effect)
dose_drug <- ev(amt = dose_nmol, cmt = "DEPOT", time = 0)
dose_inj  <- ev(amt = 1, cmt = "INJ", time = 0)
dose_combined <- dose_drug + dose_inj

# --- 2. TEST PK MODEL ---
mod_pk <- mread("test_pk", "MODEL/PK") 

sim_pk <- mod_pk %>%
  ev(dose_drug) %>%
  mrgsim(end = 28, delta = 0.1)

# plot(sim_pk, C_plasma ~ time)
p1 <- ggplot(as_tibble(sim_pk), aes(x = time, y = C_plasma)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "PK: Plasma Concentration", 
       subtitle = "For testing: assumes the physiology is a static set of volumes and flows (no trafficking/binding/t-cells/b-cells)",
       y = "Conc (nM)", x = "Time (Days)") +
  theme_bw()
print(p1)

# --- 3. TEST TRAFFICKING MODEL ---
mod_traff <- mread("test_traff", "MODEL/PD/Trafficking")

sim_traff <- mod_traff %>%
  ev(dose_combined) %>%
  init(
    # Cells
    TC_BLOOD = total_TC, TC_SPLEEN = total_TC * 30,
    BC_BLOOD = total_BC, BC_SPLEEN = total_BC * 30,
    # Tissues (Good practice to include these too)
    TC_NODE = total_TC * 10, BC_NODE = total_BC * 10,
    TC_LYMPH = total_TC * 9, BC_LYMPH = total_BC * 9,
    # Helpers: Production needs to be on from the start
    AF_TC = 1, 
    AF_BC = 1,
    INJ = 0
  ) %>%
  mrgsim(end = 28, delta = 0.1)

# plot(sim_traff, TC_BLOOD ~ time)
p_count <- ggplot(as_tibble(sim_traff), aes(x = time, y = TC_BLOOD)) +
  geom_line(color = "red", linewidth = 1) +
  labs(
    title = "Trafficking: Blood T-Cell Count",
    subtitle = "Transient reduction due to injection-induced margination",
    y = "Cells (Count)", 
    x = "Time (Days)"
  ) +
  theme_bw()

print(p_count)

if("traf_TC_blood_spleen" %in% names(sim_traff)) {
  p_rate <- ggplot(as_tibble(sim_traff), aes(x = time, y = traf_TC_blood_spleen)) +
    geom_line(color = "purple", linewidth = 1) +
    labs(
      title = "Mechanism: Blood-to-Spleen Trafficking Rate",
      subtitle = "Spike indicates the 'Injection Effect' forcing cells out",
      y = "Rate (Cells/Day)", 
      x = "Time (Days)"
    ) +
    coord_cartesian(xlim = c(0, 5)) + # Zoom in on the first 5 days
    theme_bw()
  
  print(p_rate)
}


# --- 4. TEST BINDING MODEL ---
# # --- DEBUG TEST: ZERO BINDING ---
# mod_bind <- mread("test_binding", "MODEL/PD/Binding")
# 
# # Force binding to zero
# mod_debug <- param(mod_bind, konCD3 = 0, konCD20 = 0)
# 
# sim_debug <- mod_debug %>%
#   ev(dose_combined) %>%
#   init(
#     TC_BLOOD = total_TC, BC_BLOOD = total_BC,
#     TC_SPLEEN = total_TC * 30, BC_SPLEEN = total_BC * 30,
#     TC_NODE = total_TC * 0.03, BC_NODE = total_BC * 0.03,
#     TC_LYMPH = total_TC * 19, BC_LYMPH = total_BC * 19,
#     AF_TC = 1, AF_BC = 1, INJ = 0
#   ) %>%
#   mrgsim(end = 28, delta = 0.1)
# 
# print("Simulation with kon=0 finished!")
# plot(sim_debug, C_plasma ~ time)


# # debug
# mod_bind <- mread("test_binding", "MODEL/PD/Binding")
# sim_debug <- mod_bind %>%
#   ev(ev(amt=0, cmt="DEPOT")) %>% # Zero dose
#   init(
#     TC_BLOOD = total_TC, BC_BLOOD = total_BC, 
#     TC_SPLEEN = total_TC * 30, BC_SPLEEN = total_BC * 30, 
#     TC_NODE = total_TC * 0.03, BC_NODE = total_BC * 0.03, 
#     TC_LYMPH = total_TC * 19, BC_LYMPH = total_BC * 19, 
#     AF_TC = 1, AF_BC = 1, INJ = 0
#   ) %>%
#   mrgsim(end = 28, delta = 0.1)

# actual
mod_bind <- mread("test_binding", "MODEL/PD/Binding", preclean = TRUE)
mod_bind <- param(mod_bind, knl = 1672)
sim_bind <- mod_bind %>%
  ev(dose_combined) %>%
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
  mrgsim(end = 28, delta = 0.1)

# Plot Blood vs Spleen Trimer
# plot(sim_bind, C_Trimer_BL ~ time)
# plot(sim_bind, C_Trimer_SP ~ time)
p_blood <- ggplot(as_tibble(sim_bind), aes(x = time)) +
geom_line(aes(y = C_Trimer_BL, color = "Blood"), linewidth = 1) +
scale_color_manual(name = "Tissue",
                   values = c("Blood" = "darkgreen")) +
labs(
  title = "Pharmacodynamics: Trimer Formation",
  subtitle = "Comparison of drug-target engagement in Blood",
  y = "Trimer Concentration (nM)",
  x = "Time (Days)"
) +
theme_bw() +
theme(legend.position = "right")

p_spleen <- ggplot(as_tibble(sim_bind), aes(x = time)) +
  geom_line(aes(y = C_Trimer_SP, color = "Spleen"), linewidth = 1) +
  scale_color_manual(name = "Tissue",
                     values = c("Spleen" = "purple")) +
  labs(
    title = "Pharmacodynamics: Trimer Formation",
    subtitle = "Comparison of drug-target engagement in Spleen",
    y = "Trimer Concentration (nM)",
    x = "Time (Days)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p_blood)
print(p_spleen)

# --- 5. TEST ACTIVATION MODEL ---
mod_act <- mread("test_activation", "MODEL/PD/T-Cell-Activation", preclean = TRUE)

sim_act <- mod_act %>%
  ev(dose_combined) %>%
  init(
    TC_BLOOD = total_TC, BC_BLOOD = total_BC, 
    TC_SPLEEN = total_TC * 30, BC_SPLEEN = total_BC * 30, 
    TC_NODE = total_TC * 0.03, BC_NODE = total_BC * 0.03, 
    TC_LYMPH = total_TC * 19, BC_LYMPH = total_BC * 19, 
    AF_TC = 1, AF_BC = 1, INJ = 0,
    # Activation starts at 0 (handled in $MAIN, but being explicit is good)
    pATC_BLOOD = 0
  ) %>%
  mrgsim(end = 28, delta = 0.1, 
         maxsteps = 50000) %>% # Standard settings should be fine!
  as_tibble()

# Plot: Expansion of T-Cells
ggplot(sim_act, aes(x = time)) +
  geom_line(aes(y = pATC_BLOOD, color = "Proliferating ATC (Blood)"), linewidth = 1) +
  labs(title = "T-Cell Activation", y = "Cells (Count)", x = "Time (Days)") +
  theme_bw()