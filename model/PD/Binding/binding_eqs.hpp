// BINDING EQUATIONS

// === A. BLOOD ===
double C_FreeCD3_BL   = FREE_CD3_BLOOD / Vplasma;
double C_FreeCD20_BL  = FREE_CD20_BLOOD / Vplasma;
double C_DimerCD3_BL  = DIMER_CD3_BLOOD / Vplasma;
double C_DimerCD20_BL = DIMER_CD20_BLOOD / Vplasma;
double C_Trimer_BL    = TRIMER_BLOOD / Vplasma;

double vf_DimerCD3_BL  = konCD3 * C_plasma * C_FreeCD3_BL;
double vf_DimerCD20_BL = konCD20 * C_plasma * C_FreeCD20_BL;
double vb_DimerCD3_BL  = koffCD3 * C_DimerCD3_BL;
double vb_DimerCD20_BL = koffCD20 * C_DimerCD20_BL;
double vf_Tri_viaCD3_BL  = konCD20 * C_DimerCD3_BL * C_FreeCD20_BL;
double vf_Tri_viaCD20_BL = konCD3 * C_DimerCD20_BL * C_FreeCD3_BL;
double vb_Tri_toCD3_BL   = koffCD20 * C_Trimer_BL;
double vb_Tri_toCD20_BL  = koffCD3 * C_Trimer_BL;

double syn_CD3_BL  = kdegCD3 * TC_BLOOD * R_CD3 * nmol_per_molecule;
double syn_CD20_BL = kdegCD20 * BC_BLOOD * R_CD20 * nmol_per_molecule;

dxdt_FREE_CD3_BLOOD   = -vf_DimerCD3_BL + vb_DimerCD3_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL - kdegCD3 * C_FreeCD3_BL * Vplasma + syn_CD3_BL;
dxdt_FREE_CD20_BLOOD  = -vf_DimerCD20_BL + vb_DimerCD20_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL - kdegCD20 * C_FreeCD20_BL * Vplasma + syn_CD20_BL;
dxdt_DIMER_CD3_BLOOD  = vf_DimerCD3_BL - vb_DimerCD3_BL - vf_Tri_viaCD3_BL + vb_Tri_toCD3_BL + vb_Tri_toCD20_BL - kintCD3 * C_DimerCD3_BL * Vplasma;
dxdt_DIMER_CD20_BLOOD = vf_DimerCD20_BL - vb_DimerCD20_BL - vf_Tri_viaCD20_BL + vb_Tri_toCD20_BL + vb_Tri_toCD3_BL - kintCD20 * C_DimerCD20_BL * Vplasma;
dxdt_TRIMER_BLOOD     = vf_Tri_viaCD3_BL + vf_Tri_viaCD20_BL - vb_Tri_toCD3_BL - vb_Tri_toCD20_BL;

// Update PK (Consumption)
dxdt_PLASMA = dxdt_PLASMA - (vf_DimerCD3_BL - vb_DimerCD3_BL)*Vplasma - (vf_DimerCD20_BL - vb_DimerCD20_BL)*Vplasma;

// === B. SPLEEN ===
double C_FreeCD3_SP   = FREE_CD3_SPLEEN / Vspleen;
double C_FreeCD20_SP  = FREE_CD20_SPLEEN / Vspleen;
double C_DimerCD3_SP  = DIMER_CD3_SPLEEN / Vspleen;
double C_DimerCD20_SP = DIMER_CD20_SPLEEN / Vspleen;
double C_Trimer_SP    = TRIMER_SPLEEN / Vspleen;

double vf_DimerCD3_SP  = konCD3 * C_spleen * C_FreeCD3_SP; // Uses C_spleen!
double vf_DimerCD20_SP = konCD20 * C_spleen * C_FreeCD20_SP;
double vb_DimerCD3_SP  = koffCD3 * C_DimerCD3_SP;
double vb_DimerCD20_SP = koffCD20 * C_DimerCD20_SP;
double vf_Tri_viaCD3_SP  = konCD20 * C_DimerCD3_SP * C_FreeCD20_SP;
double vf_Tri_viaCD20_SP = konCD3 * C_DimerCD20_SP * C_FreeCD3_SP;
double vb_Tri_toCD3_SP   = koffCD20 * C_Trimer_SP;
double vb_Tri_toCD20_SP  = koffCD3 * C_Trimer_SP;

double syn_CD3_SP  = kdegCD3 * TC_SPLEEN * R_CD3 * nmol_per_molecule;
double syn_CD20_SP = kdegCD20 * BC_SPLEEN * R_CD20 * nmol_per_molecule;

dxdt_FREE_CD3_SPLEEN   = -vf_DimerCD3_SP + vb_DimerCD3_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP - kdegCD3 * C_FreeCD3_SP * Vspleen + syn_CD3_SP;
dxdt_FREE_CD20_SPLEEN  = -vf_DimerCD20_SP + vb_DimerCD20_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP - kdegCD20 * C_FreeCD20_SP * Vspleen + syn_CD20_SP;
dxdt_DIMER_CD3_SPLEEN  = vf_DimerCD3_SP - vb_DimerCD3_SP - vf_Tri_viaCD3_SP + vb_Tri_toCD3_SP + vb_Tri_toCD20_SP - kintCD3 * C_DimerCD3_SP * Vspleen;
dxdt_DIMER_CD20_SPLEEN = vf_DimerCD20_SP - vb_DimerCD20_SP - vf_Tri_viaCD20_SP + vb_Tri_toCD20_SP + vb_Tri_toCD3_SP - kintCD20 * C_DimerCD20_SP * Vspleen;
dxdt_TRIMER_SPLEEN     = vf_Tri_viaCD3_SP + vf_Tri_viaCD20_SP - vb_Tri_toCD3_SP - vb_Tri_toCD20_SP;

// Update PK (Consumption in Spleen)
dxdt_SPLEEN = dxdt_SPLEEN - (vf_DimerCD3_SP - vb_DimerCD3_SP)*Vspleen - (vf_DimerCD20_SP - vb_DimerCD20_SP)*Vspleen;

// === C. LYMPH NODE ===
double C_FreeCD3_LN   = FREE_CD3_NODE / Vnode;
double C_FreeCD20_LN  = FREE_CD20_NODE / Vnode;
double C_DimerCD3_LN  = DIMER_CD3_NODE / Vnode;
double C_DimerCD20_LN = DIMER_CD20_NODE / Vnode;
double C_Trimer_LN    = TRIMER_NODE / Vnode;

double vf_DimerCD3_LN  = konCD3 * C_node * C_FreeCD3_LN; // Uses C_node!
double vf_DimerCD20_LN = konCD20 * C_node * C_FreeCD20_LN;
double vb_DimerCD3_LN  = koffCD3 * C_DimerCD3_LN;
double vb_DimerCD20_LN = koffCD20 * C_DimerCD20_LN;
double vf_Tri_viaCD3_LN  = konCD20 * C_DimerCD3_LN * C_FreeCD20_LN;
double vf_Tri_viaCD20_LN = konCD3 * C_DimerCD20_LN * C_FreeCD3_LN;
double vb_Tri_toCD3_LN   = koffCD20 * C_Trimer_LN;
double vb_Tri_toCD20_LN  = koffCD3 * C_Trimer_LN;

double syn_CD3_LN  = kdegCD3 * TC_NODE * R_CD3 * nmol_per_molecule;
double syn_CD20_LN = kdegCD20 * BC_NODE * R_CD20 * nmol_per_molecule;

dxdt_FREE_CD3_NODE   = -vf_DimerCD3_LN + vb_DimerCD3_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN - kdegCD3 * C_FreeCD3_LN * Vnode + syn_CD3_LN;
dxdt_FREE_CD20_NODE  = -vf_DimerCD20_LN + vb_DimerCD20_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN - kdegCD20 * C_FreeCD20_LN * Vnode + syn_CD20_LN;
dxdt_DIMER_CD3_NODE  = vf_DimerCD3_LN - vb_DimerCD3_LN - vf_Tri_viaCD3_LN + vb_Tri_toCD3_LN + vb_Tri_toCD20_LN - kintCD3 * C_DimerCD3_LN * Vnode;
dxdt_DIMER_CD20_NODE = vf_DimerCD20_LN - vb_DimerCD20_LN - vf_Tri_viaCD20_LN + vb_Tri_toCD20_LN + vb_Tri_toCD3_LN - kintCD20 * C_DimerCD20_LN * Vnode;
dxdt_TRIMER_NODE     = vf_Tri_viaCD3_LN + vf_Tri_viaCD20_LN - vb_Tri_toCD3_LN - vb_Tri_toCD20_LN;

dxdt_NODE = dxdt_NODE - (vf_DimerCD3_LN - vb_DimerCD3_LN)*Vnode - (vf_DimerCD20_LN - vb_DimerCD20_LN)*Vnode;

// === D. LYMPH ===
double C_FreeCD3_LY   = FREE_CD3_LYMPH / Vlymph;
double C_FreeCD20_LY  = FREE_CD20_LYMPH / Vlymph;
double C_DimerCD3_LY  = DIMER_CD3_LYMPH / Vlymph;
double C_DimerCD20_LY = DIMER_CD20_LYMPH / Vlymph;
double C_Trimer_LY    = TRIMER_LYMPH / Vlymph;

double vf_DimerCD3_LY  = konCD3 * C_lymph * C_FreeCD3_LY; // Uses C_lymph!
double vf_DimerCD20_LY = konCD20 * C_lymph * C_FreeCD20_LY;
double vb_DimerCD3_LY  = koffCD3 * C_DimerCD3_LY;
double vb_DimerCD20_LY = koffCD20 * C_DimerCD20_LY;
double vf_Tri_viaCD3_LY  = konCD20 * C_DimerCD3_LY * C_FreeCD20_LY;
double vf_Tri_viaCD20_LY = konCD3 * C_DimerCD20_LY * C_FreeCD3_LY;
double vb_Tri_toCD3_LY   = koffCD20 * C_Trimer_LY;
double vb_Tri_toCD20_LY  = koffCD3 * C_Trimer_LY;

double syn_CD3_LY  = kdegCD3 * TC_LYMPH * R_CD3 * nmol_per_molecule;
double syn_CD20_LY = kdegCD20 * BC_LYMPH * R_CD20 * nmol_per_molecule;

dxdt_FREE_CD3_LYMPH   = -vf_DimerCD3_LY + vb_DimerCD3_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY - kdegCD3 * C_FreeCD3_LY * Vlymph + syn_CD3_LY;
dxdt_FREE_CD20_LYMPH  = -vf_DimerCD20_LY + vb_DimerCD20_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY - kdegCD20 * C_FreeCD20_LY * Vlymph + syn_CD20_LY;
dxdt_DIMER_CD3_LYMPH  = vf_DimerCD3_LY - vb_DimerCD3_LY - vf_Tri_viaCD3_LY + vb_Tri_toCD3_LY + vb_Tri_toCD20_LY - kintCD3 * C_DimerCD3_LY * Vlymph;
dxdt_DIMER_CD20_LYMPH = vf_DimerCD20_LY - vb_DimerCD20_LY - vf_Tri_viaCD20_LY + vb_Tri_toCD20_LY + vb_Tri_toCD3_LY - kintCD20 * C_DimerCD20_LY * Vlymph;
dxdt_TRIMER_LYMPH     = vf_Tri_viaCD3_LY + vf_Tri_viaCD20_LY - vb_Tri_toCD3_LY - vb_Tri_toCD20_LY;

dxdt_LYMPH = dxdt_LYMPH - (vf_DimerCD3_LY - vb_DimerCD3_LY)*Vlymph - (vf_DimerCD20_LY - vb_DimerCD20_LY)*Vlymph;