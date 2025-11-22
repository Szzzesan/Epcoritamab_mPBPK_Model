// === A. BLOOD ===
dxdt_FREE_CD3_BLOOD   = -vf_DimerCD3_BL + vb_DimerCD3_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL - kdegCD3 * C_FreeCD3_BL * Vplasma + syn_CD3_BL;
dxdt_FREE_CD20_BLOOD  = -vf_DimerCD20_BL + vb_DimerCD20_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL - kdegCD20 * C_FreeCD20_BL * Vplasma + syn_CD20_BL;
dxdt_DIMER_CD3_BLOOD  = vf_DimerCD3_BL - vb_DimerCD3_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL + vb_Tri_toCD20_BL - kintCD3 * C_DimerCD3_BL * Vplasma - koutTC * C_DimerCD3_BL * Vplasma;
dxdt_DIMER_CD20_BLOOD = vf_DimerCD20_BL - vb_DimerCD20_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL + vb_Tri_toCD3_BL - kintCD20 * C_DimerCD20_BL * Vplasma - koutBC * C_DimerCD20_BL * Vplasma;
dxdt_TRIMER_BLOOD     = vf_Tri_viaCD3_BL + vf_Tri_viaCD20_BL - vb_Tri_toCD3_BL - vb_Tri_toCD20_BL - (koutTC + koutBC) * C_Trimer_BL * Vplasma;

// Update PK (Consumption)
dxdt_PLASMA = dxdt_PLASMA - (vf_DimerCD3_BL - vb_DimerCD3_BL) - (vf_DimerCD20_BL - vb_DimerCD20_BL);

// === B. SPLEEN ===
dxdt_FREE_CD3_SPLEEN   = -vf_DimerCD3_SP + vb_DimerCD3_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP - kdegCD3 * C_FreeCD3_SP * Vspleen + syn_CD3_SP;
dxdt_FREE_CD20_SPLEEN  = -vf_DimerCD20_SP + vb_DimerCD20_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP - kdegCD20 * C_FreeCD20_SP * Vspleen + syn_CD20_SP;
dxdt_DIMER_CD3_SPLEEN  = vf_DimerCD3_SP - vb_DimerCD3_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP + vb_Tri_toCD20_SP - kintCD3 * C_DimerCD3_SP * Vspleen - koutTC * C_DimerCD3_SP * Vspleen;
dxdt_DIMER_CD20_SPLEEN = vf_DimerCD20_SP - vb_DimerCD20_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP + vb_Tri_toCD3_SP - kintCD20 * C_DimerCD20_SP * Vspleen - koutBC * C_DimerCD20_SP * Vspleen;
dxdt_TRIMER_SPLEEN     = vf_Tri_viaCD3_SP + vf_Tri_viaCD20_SP - vb_Tri_toCD3_SP - vb_Tri_toCD20_SP - (koutTC + koutBC) * C_Trimer_SP * Vspleen;

// Update PK (Consumption in Spleen)
dxdt_SPLEEN = dxdt_SPLEEN - (vf_DimerCD3_SP - vb_DimerCD3_SP) - (vf_DimerCD20_SP - vb_DimerCD20_SP);

// === C. LYMPH NODE ===
dxdt_FREE_CD3_NODE   = -vf_DimerCD3_LN + vb_DimerCD3_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN - kdegCD3 * C_FreeCD3_LN * Vnode + syn_CD3_LN;
dxdt_FREE_CD20_NODE  = -vf_DimerCD20_LN + vb_DimerCD20_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN - kdegCD20 * C_FreeCD20_LN * Vnode + syn_CD20_LN;
dxdt_DIMER_CD3_NODE  = vf_DimerCD3_LN - vb_DimerCD3_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN + vb_Tri_toCD20_LN - kintCD3 * C_DimerCD3_LN * Vnode - koutTC * C_DimerCD3_LN * Vnode;
dxdt_DIMER_CD20_NODE = vf_DimerCD20_LN - vb_DimerCD20_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN + vb_Tri_toCD3_LN - kintCD20 * C_DimerCD20_LN * Vnode - koutBC * C_DimerCD20_LN * Vnode;
dxdt_TRIMER_NODE     = vf_Tri_viaCD3_LN + vf_Tri_viaCD20_LN - vb_Tri_toCD3_LN - vb_Tri_toCD20_LN - (koutTC + koutBC) * C_Trimer_LN * Vnode;

dxdt_NODE = dxdt_NODE - (vf_DimerCD3_LN - vb_DimerCD3_LN) - (vf_DimerCD20_LN - vb_DimerCD20_LN);

// === D. LYMPH ===
dxdt_FREE_CD3_LYMPH   = -vf_DimerCD3_LY + vb_DimerCD3_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY - kdegCD3 * C_FreeCD3_LY * Vlymph + syn_CD3_LY;
dxdt_FREE_CD20_LYMPH  = -vf_DimerCD20_LY + vb_DimerCD20_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY - kdegCD20 * C_FreeCD20_LY * Vlymph + syn_CD20_LY;
dxdt_DIMER_CD3_LYMPH  = vf_DimerCD3_LY - vb_DimerCD3_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY + vb_Tri_toCD20_LY - kintCD3 * C_DimerCD3_LY * Vlymph - koutTC * C_DimerCD3_LY * Vlymph;
dxdt_DIMER_CD20_LYMPH = vf_DimerCD20_LY - vb_DimerCD20_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY + vb_Tri_toCD3_LY - kintCD20 * C_DimerCD20_LY * Vlymph - koutBC * C_DimerCD20_LY * Vlymph;
dxdt_TRIMER_LYMPH     = vf_Tri_viaCD3_LY + vf_Tri_viaCD20_LY - vb_Tri_toCD3_LY - vb_Tri_toCD20_LY - (koutTC + koutBC) * C_Trimer_LY * Vlymph;

dxdt_LYMPH = dxdt_LYMPH - (vf_DimerCD3_LY - vb_DimerCD3_LY) - (vf_DimerCD20_LY - vb_DimerCD20_LY);