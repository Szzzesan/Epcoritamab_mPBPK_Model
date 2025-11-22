// --- vATC_BC (Traffics like T-Cells) ---
dxdt_vATC_BC_BLOOD  = rate_act_BC_BL + flow_vATC_LY_BL - flow_vATC_BL_SP - rate_death_ATC * vATC_BC_BLOOD;
dxdt_vATC_BC_SPLEEN = rate_act_BC_SP + flow_vATC_BL_SP - flow_vATC_SP_LN - rate_death_ATC * vATC_BC_SPLEEN;
dxdt_vATC_BC_NODE   = rate_act_BC_LN + flow_vATC_SP_LN - flow_vATC_LN_LY - rate_death_ATC * vATC_BC_NODE;
dxdt_vATC_BC_LYMPH  = rate_act_BC_LY + flow_vATC_LN_LY - flow_vATC_LY_BL - rate_death_ATC * vATC_BC_LYMPH;

// --- vATC_TUMOR (Immobile - Node Only) ---
dxdt_vATC_TUMOR_NODE = rate_act_Tumor - rate_death_ATC * vATC_TUMOR_NODE;

// --- pATC (Traffics like T-Cells) ---
dxdt_pATC_BLOOD  = exp_BL_rate + flow_pATC_LY_BL - flow_pATC_BL_SP - rate_death_ATC * pATC_BLOOD;
dxdt_pATC_SPLEEN = exp_SP_rate + flow_pATC_BL_SP - flow_pATC_SP_LN - rate_death_ATC * pATC_SPLEEN;
dxdt_pATC_NODE   = exp_LN_rate + flow_pATC_SP_LN - flow_pATC_LN_LY - rate_death_ATC * pATC_NODE;
dxdt_pATC_LYMPH  = exp_LY_rate + flow_pATC_LN_LY - flow_pATC_LY_BL - rate_death_ATC * pATC_LYMPH;