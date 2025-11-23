// Units should be nmol/day
// Global Scalar to handle stiffness
// This slows reaction time from ~seconds to ~minutes, helping lsoda
double scale_binding = 1.0;

// === A. BLOOD ===
double C_FreeCD3_BL   = FREE_CD3_BLOOD / Vplasma;
double C_FreeCD20_BL  = FREE_CD20_BLOOD / Vplasma;
double C_DimerCD3_BL  = DIMER_CD3_BLOOD / Vplasma;
double C_DimerCD20_BL = DIMER_CD20_BLOOD / Vplasma;
double C_Trimer_BL    = TRIMER_BLOOD / Vplasma;
// double C_FreeCD3_BL   = CLAMP(FREE_CD3_BLOOD) / Vplasma;
// double C_FreeCD20_BL  = CLAMP(FREE_CD20_BLOOD) / Vplasma;
// double C_DimerCD3_BL  = CLAMP(DIMER_CD3_BLOOD) / Vplasma;
// double C_DimerCD20_BL = CLAMP(DIMER_CD20_BLOOD) / Vplasma;
// double C_Trimer_BL    = CLAMP(TRIMER_BLOOD) / Vplasma;

double vf_DimerCD3_BL  = scale_binding * konCD3 * C_plasma * C_FreeCD3_BL * Vplasma;
double vf_DimerCD20_BL = scale_binding * konCD20 * C_plasma * C_FreeCD20_BL * Vplasma;
double vb_DimerCD3_BL  = scale_binding * koffCD3 * C_DimerCD3_BL * Vplasma;
double vb_DimerCD20_BL = scale_binding * koffCD20 * C_DimerCD20_BL * Vplasma;
double vf_Tri_viaCD3_BL  = scale_binding * konCD20 * C_DimerCD3_BL * C_FreeCD20_BL * Vplasma;
double vf_Tri_viaCD20_BL = scale_binding * konCD3 * C_DimerCD20_BL * C_FreeCD3_BL * Vplasma;
double vb_Tri_toCD3_BL   = scale_binding * koffCD20 * C_Trimer_BL * Vplasma;
double vb_Tri_toCD20_BL  = scale_binding * koffCD3 * C_Trimer_BL * Vplasma;

double syn_CD3_BL  = kdegCD3 * TC_BLOOD * R_CD3 * nmol_per_molecule;
double syn_CD20_BL = kdegCD20 * BC_BLOOD * R_CD20 * nmol_per_molecule;

// === B. SPLEEN ===
double C_FreeCD3_SP   = FREE_CD3_SPLEEN / Vspleen;
double C_FreeCD20_SP  = FREE_CD20_SPLEEN / Vspleen;
double C_DimerCD3_SP  = DIMER_CD3_SPLEEN / Vspleen;
double C_DimerCD20_SP = DIMER_CD20_SPLEEN / Vspleen;
double C_Trimer_SP    = TRIMER_SPLEEN / Vspleen;
// double C_FreeCD3_SP   = CLAMP(FREE_CD3_SPLEEN) / Vspleen;
// double C_FreeCD20_SP  = CLAMP(FREE_CD20_SPLEEN) / Vspleen;
// double C_DimerCD3_SP  = CLAMP(DIMER_CD3_SPLEEN) / Vspleen;
// double C_DimerCD20_SP = CLAMP(DIMER_CD20_SPLEEN) / Vspleen;
// double C_Trimer_SP    = CLAMP(TRIMER_SPLEEN) / Vspleen;

double vf_DimerCD3_SP  = scale_binding * konCD3 * C_spleen * C_FreeCD3_SP * Vspleen;
double vf_DimerCD20_SP = scale_binding * konCD20 * C_spleen * C_FreeCD20_SP * Vspleen;
double vb_DimerCD3_SP  = scale_binding * koffCD3 * C_DimerCD3_SP * Vspleen;
double vb_DimerCD20_SP = scale_binding * koffCD20 * C_DimerCD20_SP * Vspleen;
double vf_Tri_viaCD3_SP  = scale_binding * konCD20 * C_DimerCD3_SP * C_FreeCD20_SP * Vspleen;
double vf_Tri_viaCD20_SP = scale_binding * konCD3 * C_DimerCD20_SP * C_FreeCD3_SP * Vspleen;
double vb_Tri_toCD3_SP   = scale_binding * koffCD20 * C_Trimer_SP * Vspleen;
double vb_Tri_toCD20_SP  = scale_binding * koffCD3 * C_Trimer_SP * Vspleen;

double syn_CD3_SP  = kdegCD3 * TC_SPLEEN * R_CD3 * nmol_per_molecule;
double syn_CD20_SP = kdegCD20 * BC_SPLEEN * R_CD20 * nmol_per_molecule;

// === C. LYMPH NODE ===
double C_FreeCD3_LN   = FREE_CD3_NODE / Vnode;
double C_FreeCD20_LN  = FREE_CD20_NODE / Vnode;
double C_DimerCD3_LN  = DIMER_CD3_NODE / Vnode;
double C_DimerCD20_LN = DIMER_CD20_NODE / Vnode;
double C_Trimer_LN    = TRIMER_NODE / Vnode;
// double C_FreeCD3_LN   = CLAMP(FREE_CD3_NODE) / Vnode;
// double C_FreeCD20_LN  = CLAMP(FREE_CD20_NODE) / Vnode;
// double C_DimerCD3_LN  = CLAMP(DIMER_CD3_NODE) / Vnode;
// double C_DimerCD20_LN = CLAMP(DIMER_CD20_NODE) / Vnode;
// double C_Trimer_LN    = CLAMP(TRIMER_NODE) / Vnode;

double vf_DimerCD3_LN  = scale_binding * konCD3 * C_node * C_FreeCD3_LN * Vnode;
double vf_DimerCD20_LN = scale_binding * konCD20 * C_node * C_FreeCD20_LN * Vnode;
double vb_DimerCD3_LN  = scale_binding * koffCD3 * C_DimerCD3_LN * Vnode;
double vb_DimerCD20_LN = scale_binding * koffCD20 * C_DimerCD20_LN * Vnode;
double vf_Tri_viaCD3_LN  = scale_binding * konCD20 * C_DimerCD3_LN * C_FreeCD20_LN * Vnode;
double vf_Tri_viaCD20_LN = scale_binding * konCD3 * C_DimerCD20_LN * C_FreeCD3_LN * Vnode;
double vb_Tri_toCD3_LN   = scale_binding * koffCD20 * C_Trimer_LN * Vnode;
double vb_Tri_toCD20_LN  = scale_binding * koffCD3 * C_Trimer_LN * Vnode;

double syn_CD3_LN  = kdegCD3 * TC_NODE * R_CD3 * nmol_per_molecule;
double syn_CD20_LN = kdegCD20 * BC_NODE * R_CD20 * nmol_per_molecule;

// === D. LYMPH ===
double C_FreeCD3_LY   = FREE_CD3_LYMPH / Vlymph;
double C_FreeCD20_LY  = FREE_CD20_LYMPH / Vlymph;
double C_DimerCD3_LY  = DIMER_CD3_LYMPH / Vlymph;
double C_DimerCD20_LY = DIMER_CD20_LYMPH / Vlymph;
double C_Trimer_LY    = TRIMER_LYMPH / Vlymph;
// double C_FreeCD3_LY   = CLAMP(FREE_CD3_LYMPH) / Vlymph;
// double C_FreeCD20_LY  = CLAMP(FREE_CD20_LYMPH) / Vlymph;
// double C_DimerCD3_LY  = CLAMP(DIMER_CD3_LYMPH) / Vlymph;
// double C_DimerCD20_LY = CLAMP(DIMER_CD20_LYMPH) / Vlymph;
// double C_Trimer_LY    = CLAMP(TRIMER_LYMPH) / Vlymph;

double vf_DimerCD3_LY  = scale_binding * konCD3 * C_lymph * C_FreeCD3_LY * Vlymph;
double vf_DimerCD20_LY = scale_binding * konCD20 * C_lymph * C_FreeCD20_LY * Vlymph;
double vb_DimerCD3_LY  = scale_binding * koffCD3 * C_DimerCD3_LY * Vlymph;
double vb_DimerCD20_LY = scale_binding * koffCD20 * C_DimerCD20_LY * Vlymph;
double vf_Tri_viaCD3_LY  = scale_binding * konCD20 * C_DimerCD3_LY * C_FreeCD20_LY * Vlymph;
double vf_Tri_viaCD20_LY = scale_binding * konCD3 * C_DimerCD20_LY * C_FreeCD3_LY * Vlymph;
double vb_Tri_toCD3_LY   = scale_binding * koffCD20 * C_Trimer_LY * Vlymph;
double vb_Tri_toCD20_LY  = scale_binding * koffCD3 * C_Trimer_LY * Vlymph;

double syn_CD3_LY  = kdegCD3 * TC_LYMPH * R_CD3 * nmol_per_molecule;
double syn_CD20_LY = kdegCD20 * BC_LYMPH * R_CD20 * nmol_per_molecule;