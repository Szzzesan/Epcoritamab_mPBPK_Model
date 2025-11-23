// --- TRAFFICKING RATES  ---
// Rate Constants
double rate_BL_SP = kpt * (1 + INJ_Scaler * INJ); // Blood to Spleen (with Injection Effect)
double rate_SP_LN = ktn;                          // Spleen to Node
double rate_LN_LY = knl;                           // knl
double rate_LY_BL = klp;                          // Lymph to Blood

// === A. BLOOD ===
double flow_FreeCD3_BL   = rate_LY_BL * FREE_CD3_LYMPH   - rate_BL_SP * FREE_CD3_BLOOD;
double flow_FreeCD20_BL  = rate_LY_BL * FREE_CD20_LYMPH - rate_BL_SP * FREE_CD20_BLOOD;
double flow_DimerCD3_BL  = rate_LY_BL * DIMER_CD3_LYMPH - rate_BL_SP * DIMER_CD3_BLOOD;
double flow_DimerCD20_BL = rate_LY_BL * DIMER_CD20_LYMPH- rate_BL_SP * DIMER_CD20_BLOOD;
double flow_Trimer_BL    = rate_LY_BL * TRIMER_LYMPH    - rate_BL_SP * TRIMER_BLOOD;

// 1. Free Receptors
dxdt_FREE_CD3_BLOOD   = -vf_DimerCD3_BL + vb_DimerCD3_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL 
                        - kdegCD3 * C_FreeCD3_BL * Vplasma + syn_CD3_BL + flow_FreeCD3_BL 
                        - koutTC * C_FreeCD3_BL * Vplasma 
                        + prod_new_CD3_BL; 

dxdt_FREE_CD20_BLOOD  = -vf_DimerCD20_BL + vb_DimerCD20_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL 
                        - kdegCD20 * C_FreeCD20_BL * Vplasma + syn_CD20_BL + flow_FreeCD20_BL 
                        - (koutBC + rate_kill_BC_BL) * C_FreeCD20_BL * Vplasma 
                        + prod_new_CD20_BL;
                        
dxdt_DIMER_CD3_BLOOD  = vf_DimerCD3_BL - vb_DimerCD3_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL  
                        - kintCD3 * C_DimerCD3_BL * Vplasma - koutTC * C_DimerCD3_BL * Vplasma + flow_DimerCD3_BL
                        + (koutBC + rate_kill_BC_BL) * C_Trimer_BL * Vplasma;

dxdt_DIMER_CD20_BLOOD = vf_DimerCD20_BL - vb_DimerCD20_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL  
                        - kintCD20 * C_DimerCD20_BL * Vplasma - koutBC * C_DimerCD20_BL * Vplasma + flow_DimerCD20_BL
                        + koutTC * C_Trimer_BL * Vplasma - rate_kill_BC_BL * C_DimerCD20_BL * Vplasma;

dxdt_TRIMER_BLOOD     = vf_Tri_viaCD3_BL + vf_Tri_viaCD20_BL - vb_Tri_toCD3_BL - vb_Tri_toCD20_BL 
                        - (koutTC + koutBC) * C_Trimer_BL * Vplasma + flow_Trimer_BL
                        - rate_kill_BC_BL * C_Trimer_BL * Vplasma;

// Update PK (Consumption)
dxdt_PLASMA = dxdt_PLASMA - (vf_DimerCD3_BL - vb_DimerCD3_BL) - (vf_DimerCD20_BL - vb_DimerCD20_BL);

// === B. SPLEEN ===
double flow_FreeCD3_SP   = rate_BL_SP * FREE_CD3_BLOOD  - rate_SP_LN * FREE_CD3_SPLEEN;
double flow_FreeCD20_SP  = rate_BL_SP * FREE_CD20_BLOOD - rate_SP_LN * FREE_CD20_SPLEEN;
double flow_DimerCD3_SP  = rate_BL_SP * DIMER_CD3_BLOOD - rate_SP_LN * DIMER_CD3_SPLEEN;
double flow_DimerCD20_SP = rate_BL_SP * DIMER_CD20_BLOOD- rate_SP_LN * DIMER_CD20_SPLEEN;
double flow_Trimer_SP    = rate_BL_SP * TRIMER_BLOOD    - rate_SP_LN * TRIMER_SPLEEN;

dxdt_FREE_CD3_SPLEEN   = -vf_DimerCD3_SP + vb_DimerCD3_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP 
                         - kdegCD3 * C_FreeCD3_SP * Vspleen + syn_CD3_SP + flow_FreeCD3_SP
                         - koutTC * C_FreeCD3_SP * Vspleen;

dxdt_FREE_CD20_SPLEEN  = -vf_DimerCD20_SP + vb_DimerCD20_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP 
                         - kdegCD20 * C_FreeCD20_SP * Vspleen + syn_CD20_SP + flow_FreeCD20_SP
                         - (koutBC + rate_kill_BC_SP) * C_FreeCD20_SP * Vspleen;

dxdt_DIMER_CD3_SPLEEN  = vf_DimerCD3_SP - vb_DimerCD3_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP + vb_Tri_toCD20_SP 
                         - kintCD3 * C_DimerCD3_SP * Vspleen - koutTC * C_DimerCD3_SP * Vspleen + flow_DimerCD3_SP
                         + (koutBC + rate_kill_BC_SP) * C_Trimer_SP * Vspleen;

dxdt_DIMER_CD20_SPLEEN = vf_DimerCD20_SP - vb_DimerCD20_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP + vb_Tri_toCD3_SP 
                         - kintCD20 * C_DimerCD20_SP * Vspleen - koutBC * C_DimerCD20_SP * Vspleen + flow_DimerCD20_SP
                         + koutTC * C_Trimer_SP * Vspleen - rate_kill_BC_SP * C_DimerCD20_SP * Vspleen;

dxdt_TRIMER_SPLEEN     = vf_Tri_viaCD3_SP + vf_Tri_viaCD20_SP - vb_Tri_toCD3_SP - vb_Tri_toCD20_SP 
                         - (koutTC + koutBC) * C_Trimer_SP * Vspleen + flow_Trimer_SP
                         - rate_kill_BC_SP * C_Trimer_SP * Vspleen;

// Update PK (Consumption in Spleen)
dxdt_SPLEEN = dxdt_SPLEEN - (vf_DimerCD3_SP - vb_DimerCD3_SP) - (vf_DimerCD20_SP - vb_DimerCD20_SP);

// === C. LYMPH NODE ===
double flow_FreeCD3_LN   = rate_SP_LN * FREE_CD3_SPLEEN - rate_LN_LY * FREE_CD3_NODE;
double flow_FreeCD20_LN  = rate_SP_LN * FREE_CD20_SPLEEN- rate_LN_LY * FREE_CD20_NODE;
double flow_DimerCD3_LN  = rate_SP_LN * DIMER_CD3_SPLEEN- rate_LN_LY * DIMER_CD3_NODE;
double flow_DimerCD20_LN = rate_SP_LN * DIMER_CD20_SPLEEN-rate_LN_LY * DIMER_CD20_NODE;
double flow_Trimer_LN    = rate_SP_LN * TRIMER_SPLEEN   - rate_LN_LY * TRIMER_NODE;

dxdt_FREE_CD3_NODE   = -vf_DimerCD3_LN + vb_DimerCD3_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN 
                       - kdegCD3 * C_FreeCD3_LN * Vnode + syn_CD3_LN + flow_FreeCD3_LN
                       - koutTC * C_FreeCD3_LN * Vnode;

dxdt_FREE_CD20_NODE  = -vf_DimerCD20_LN + vb_DimerCD20_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN 
                       - kdegCD20 * C_FreeCD20_LN * Vnode + syn_CD20_LN + flow_FreeCD20_LN
                       - (koutBC + rate_kill_BC_LN) * C_FreeCD20_LN * Vnode;

dxdt_DIMER_CD3_NODE  = vf_DimerCD3_LN - vb_DimerCD3_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN + vb_Tri_toCD20_LN 
                       - kintCD3 * C_DimerCD3_LN * Vnode - koutTC * C_DimerCD3_LN * Vnode + flow_DimerCD3_LN
                       + (koutBC + rate_kill_BC_LN) * C_Trimer_LN * Vnode;

dxdt_DIMER_CD20_NODE = vf_DimerCD20_LN - vb_DimerCD20_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN + vb_Tri_toCD3_LN 
                       - kintCD20 * C_DimerCD20_LN * Vnode - koutBC * C_DimerCD20_LN * Vnode + flow_DimerCD20_LN
                       + koutTC * C_Trimer_LN * Vnode - rate_kill_BC_LN * C_DimerCD20_LN * Vnode;

dxdt_TRIMER_NODE     = vf_Tri_viaCD3_LN + vf_Tri_viaCD20_LN - vb_Tri_toCD3_LN - vb_Tri_toCD20_LN 
                       - (koutTC + koutBC) * C_Trimer_LN * Vnode + flow_Trimer_LN
                       - rate_kill_BC_LN * C_Trimer_LN * Vnode;

dxdt_NODE = dxdt_NODE - (vf_DimerCD3_LN - vb_DimerCD3_LN) - (vf_DimerCD20_LN - vb_DimerCD20_LN);

// === D. LYMPH ===
double flow_FreeCD3_LY   = rate_LN_LY * FREE_CD3_NODE - rate_LY_BL * FREE_CD3_LYMPH;
double flow_FreeCD20_LY  = rate_LN_LY * FREE_CD20_NODE- rate_LY_BL * FREE_CD20_LYMPH;
double flow_DimerCD3_LY  = rate_LN_LY * DIMER_CD3_NODE - rate_LY_BL * DIMER_CD3_LYMPH;
double flow_DimerCD20_LY = rate_LN_LY * DIMER_CD20_NODE- rate_LY_BL * DIMER_CD20_LYMPH;
double flow_Trimer_LY    = rate_LN_LY * TRIMER_NODE    - rate_LY_BL * TRIMER_LYMPH;

dxdt_FREE_CD3_LYMPH   = -vf_DimerCD3_LY + vb_DimerCD3_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY 
                        - kdegCD3 * C_FreeCD3_LY * Vlymph + syn_CD3_LY + flow_FreeCD3_LY
                        - koutTC * C_FreeCD3_LY * Vlymph;

dxdt_FREE_CD20_LYMPH  = -vf_DimerCD20_LY + vb_DimerCD20_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY 
                        - kdegCD20 * C_FreeCD20_LY * Vlymph + syn_CD20_LY + flow_FreeCD20_LY
                        - (koutBC + rate_kill_BC_LY) * C_FreeCD20_LY * Vlymph;

dxdt_DIMER_CD3_LYMPH  = vf_DimerCD3_LY - vb_DimerCD3_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY + vb_Tri_toCD20_LY 
                        - kintCD3 * C_DimerCD3_LY * Vlymph - koutTC * C_DimerCD3_LY * Vlymph + flow_DimerCD3_LY
                        + (koutBC + rate_kill_BC_LY) * C_Trimer_LY * Vlymph;

dxdt_DIMER_CD20_LYMPH = vf_DimerCD20_LY - vb_DimerCD20_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY + vb_Tri_toCD3_LY 
                        - kintCD20 * C_DimerCD20_LY * Vlymph - koutBC * C_DimerCD20_LY * Vlymph + flow_DimerCD20_LY
                        + koutTC * C_Trimer_LY * Vlymph - rate_kill_BC_LY * C_DimerCD20_LY * Vlymph;

dxdt_TRIMER_LYMPH     = vf_Tri_viaCD3_LY + vf_Tri_viaCD20_LY - vb_Tri_toCD3_LY - vb_Tri_toCD20_LY 
                        - (koutTC + koutBC) * C_Trimer_LY * Vlymph + flow_Trimer_LY
                        - rate_kill_BC_LY * C_Trimer_LY * Vlymph;

dxdt_LYMPH = dxdt_LYMPH - (vf_DimerCD3_LY - vb_DimerCD3_LY) - (vf_DimerCD20_LY - vb_DimerCD20_LY);