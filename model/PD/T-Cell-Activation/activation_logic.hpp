// activation_logic.hpp
// This file contains a C++ function to calculate rates.
// It avoids redefinition errors by keeping variables LOCAL to the function.

#ifndef ACTIVATION_LOGIC_HPP
#define ACTIVATION_LOGIC_HPP

// Define a struct to hold all the calculated rates
struct ActivationRates {
  double d_vATC_BC_BL, d_vATC_BC_SP, d_vATC_BC_LN, d_vATC_BC_LY;
  double d_vATC_Tumor;
  double d_pATC_BL, d_pATC_SP, d_pATC_LN, d_pATC_LY;
  double total_pATC; // For output
};

// The calculation function
inline ActivationRates calculate_activation_rates(
  double TimeAfterDose,
  double BC_BLOOD, double BC_SPLEEN, double BC_NODE, double BC_LYMPH,
  double TRIMER_BLOOD, double TRIMER_SPLEEN, double TRIMER_NODE, double TRIMER_LYMPH,
  double vATC_BC_BLOOD, double vATC_BC_SPLEEN, double vATC_BC_NODE, double vATC_BC_LYMPH,
  double vATC_TUMOR_NODE,
  double pATC_BLOOD, double pATC_SPLEEN, double pATC_NODE, double pATC_LYMPH,
  double INJ,
  // Parameters
  double sim_slope, double sim_slopetumor, double expand_factor, double koutATC,
  double TAD, double Tp, double Trimer_Threshold,
  double kpt, double ktn, double knl, double klp, double INJ_Scaler,
  double nmol_per_molecule
) {
  
  // --- 1. Tumor Activation ---
  double TUMOR_CELLS = 0.0; 
  double Trimer_per_Tumor = 0.0;
  if (TUMOR_CELLS > 1.0) {
    Trimer_per_Tumor = (TRIMER_NODE > 0 ? TRIMER_NODE : 0.0) / nmol_per_molecule / TUMOR_CELLS;
  }
  
  double RELU = 0.01; 
  if (Trimer_per_Tumor > Trimer_Threshold) RELU = 1.0; 
  
  double rate_act_Tumor = 0.0;
  if (TimeAfterDose > TAD) {
    rate_act_Tumor = RELU * sim_slopetumor * Trimer_per_Tumor;
  }
  
  // --- 2. B-Cell Activation ---
  // Helper lambda for safe division
  auto get_trimer_per_bc = [&](double bc, double trimer) {
    return (bc > 1.0) ? ((trimer > 0 ? trimer : 0.0) / nmol_per_molecule) / bc : 0.0;
  };
  
  double rate_act_BC_BL = (TimeAfterDose > TAD) ? sim_slope * get_trimer_per_bc(BC_BLOOD, TRIMER_BLOOD) : 0.0;
  double rate_act_BC_SP = (TimeAfterDose > TAD) ? sim_slope * get_trimer_per_bc(BC_SPLEEN, TRIMER_SPLEEN) : 0.0;
  double rate_act_BC_LN = (TimeAfterDose > TAD) ? sim_slope * get_trimer_per_bc(BC_NODE, TRIMER_NODE) : 0.0;
  double rate_act_BC_LY = (TimeAfterDose > TAD) ? sim_slope * get_trimer_per_bc(BC_LYMPH, TRIMER_LYMPH) : 0.0;
  
  // --- 3. Trafficking Flows ---
  // Helper for clamping
  auto clamp = [](double x) { return x > 0 ? x : 0.0; };
  
  // vATC Flows
  double flow_vATC_BL_SP = kpt * (1 + INJ_Scaler * INJ) * clamp(vATC_BC_BLOOD);
  double flow_vATC_SP_LN = ktn * clamp(vATC_BC_SPLEEN);
  double flow_vATC_LN_LY = knl * clamp(vATC_BC_NODE);
  double flow_vATC_LY_BL = klp * clamp(vATC_BC_LYMPH);
  
  // pATC Flows
  double flow_pATC_BL_SP = kpt * (1 + INJ_Scaler * INJ) * clamp(pATC_BLOOD);
  double flow_pATC_SP_LN = ktn * clamp(pATC_SPLEEN);
  double flow_pATC_LN_LY = knl * clamp(pATC_NODE);
  double flow_pATC_LY_BL = klp * clamp(pATC_LYMPH);
  
  // --- 4. Expansion & Death ---
  double exp_BL_rate = expand_factor * clamp(vATC_BC_BLOOD);
  double exp_SP_rate = expand_factor * clamp(vATC_BC_SPLEEN);
  double exp_LN_rate = expand_factor * (clamp(vATC_BC_NODE) + clamp(vATC_TUMOR_NODE));
  double exp_LY_rate = expand_factor * clamp(vATC_BC_LYMPH);
  
  double rate_death_ATC = (TimeAfterDose > (TAD + Tp)) ? koutATC : 0.0;
  
  // --- 5. Calculate Derivatives ---
  ActivationRates res;
  
  res.d_vATC_BC_BL  = rate_act_BC_BL + flow_vATC_LY_BL - flow_vATC_BL_SP - rate_death_ATC * vATC_BC_BLOOD;
  res.d_vATC_BC_SP  = rate_act_BC_SP + flow_vATC_BL_SP - flow_vATC_SP_LN - rate_death_ATC * vATC_BC_SPLEEN;
  res.d_vATC_BC_LN  = rate_act_BC_LN + flow_vATC_SP_LN - flow_vATC_LN_LY - rate_death_ATC * vATC_BC_NODE;
  res.d_vATC_BC_LY  = rate_act_BC_LY + flow_vATC_LN_LY - flow_vATC_LY_BL - rate_death_ATC * vATC_BC_LYMPH;
  
  res.d_vATC_Tumor  = rate_act_Tumor - rate_death_ATC * vATC_TUMOR_NODE;
  
  res.d_pATC_BL     = exp_BL_rate + flow_pATC_LY_BL - flow_pATC_BL_SP - rate_death_ATC * pATC_BLOOD;
  res.d_pATC_SP     = exp_SP_rate + flow_pATC_BL_SP - flow_pATC_SP_LN - rate_death_ATC * pATC_SPLEEN;
  res.d_pATC_LN     = exp_LN_rate + flow_pATC_SP_LN - flow_pATC_LN_LY - rate_death_ATC * pATC_NODE;
  res.d_pATC_LY     = exp_LY_rate + flow_pATC_LN_LY - flow_pATC_LY_BL - rate_death_ATC * pATC_LYMPH;
  
  res.total_pATC    = clamp(pATC_BLOOD) + clamp(pATC_SPLEEN) + clamp(pATC_NODE) + clamp(pATC_LYMPH);
  
  return res;
}

#endif