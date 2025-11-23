// 1. Activation against Tumor (Local in Node)
// Placeholder: We haven't built the Tumor compartment yet, so assume 0 for now to avoid crash
double TUMOR_CELLS = 0.0; 
// In the next step (Killing Model), this will become a real variable.

double Trimer_per_Tumor = 0.0;

if (TUMOR_CELLS > 1.0) {
  // Logic: Trimers in Node / Number of Tumor Cells
  Trimer_per_Tumor = (TRIMER_NODE / nmol_per_molecule) / TUMOR_CELLS;
}

// Stepwise Activation Function (RELU)
double RELU = 0.01; // Baseline low activation
if (Trimer_per_Tumor > Trimer_Threshold) {
  RELU = 1.0; // Full activation if threshold met
}

double rate_act_Tumor = 0.0;
if (TimeAfterDose > TAD) {
  rate_act_Tumor = RELU * sim_slopetumor * Trimer_per_Tumor;
}

// 2. Activation against B-Cells (Systemic)
// Rate = sim_slope * (Trimers / B-Cells) in that specific tissue

double Trimer_per_BC_BL = (BC_BLOOD > 1)  ? (TRIMER_BLOOD / nmol_per_molecule) / BC_BLOOD  : 0.0;
double Trimer_per_BC_SP = (BC_SPLEEN > 1) ? (TRIMER_SPLEEN / nmol_per_molecule) / BC_SPLEEN : 0.0;
double Trimer_per_BC_LN = (BC_NODE > 1)   ? (TRIMER_NODE / nmol_per_molecule) / BC_NODE   : 0.0;
double Trimer_per_BC_LY = (BC_LYMPH > 1)  ? (TRIMER_LYMPH / nmol_per_molecule) / BC_LYMPH  : 0.0;

double rate_act_BC_BL = (TimeAfterDose > TAD) ? sim_slope * Trimer_per_BC_BL : 0.0;
double rate_act_BC_SP = (TimeAfterDose > TAD) ? sim_slope * Trimer_per_BC_SP : 0.0;
double rate_act_BC_LN = (TimeAfterDose > TAD) ? sim_slope * Trimer_per_BC_LN : 0.0;
double rate_act_BC_LY = (TimeAfterDose > TAD) ? sim_slope * Trimer_per_BC_LY : 0.0;

// 3. Trafficking in and out
// (1) vATC_BC Trafficking Flows
double flow_vATC_BL_SP = kpt * (1 + INJ_Scaler * INJ) * vATC_BC_BLOOD;
double flow_vATC_SP_LN = ktn * vATC_BC_SPLEEN;
double flow_vATC_LN_LY = knl * vATC_BC_NODE; // Use 1672
double flow_vATC_LY_BL = klp * vATC_BC_LYMPH;

// (2)pATC Flows
double flow_pATC_BL_SP = kpt * (1 + INJ_Scaler * INJ) * pATC_BLOOD;
double flow_pATC_SP_LN = ktn * pATC_SPLEEN;
double flow_pATC_LN_LY = knl * pATC_NODE;
double flow_pATC_LY_BL = klp * pATC_LYMPH;

// 4. Expansion & Death Rates
// Activated cells multiply!
// Expansion depends on local activated cells
double exp_BL_rate = expand_factor * vATC_BC_BLOOD;
double exp_SP_rate = expand_factor * vATC_BC_SPLEEN;
// Node has BOTH types of vATC contributing to expansion!
double exp_LN_rate = expand_factor * (vATC_BC_NODE + vATC_TUMOR_NODE); 
double exp_LY_rate = expand_factor * vATC_BC_LYMPH;

double rate_death_ATC = (TimeAfterDose > (TAD + Tp)) ? koutATC : 0.0;