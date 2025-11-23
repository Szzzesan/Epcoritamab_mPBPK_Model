// // T-Cells Production & Death Rates
// double prod_TC  = kinTC * pow(AF_TC, r);
// double death_TC_blood  = koutTC * TC_BLOOD;
// double death_TC_spleen = koutTC * TC_SPLEEN;
// double death_TC_node   = koutTC * TC_NODE;
// double death_TC_lymph  = koutTC * TC_LYMPH;
// 
// // T-Cells Trafficking Fluxes (Movement)
// double traf_TC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * TC_BLOOD;
// double traf_TC_spleen_node  = ktn * TC_SPLEEN;
// double traf_TC_node_lymph   = knl * TC_NODE;
// double traf_TC_lymph_blood  = klp * TC_LYMPH;
// 
// // B-Cells Production & Death Rates
// double prod_BC  = kinBC * pow(AF_BC, r);
// double death_BC_blood  = koutBC * BC_BLOOD;
// double death_BC_spleen = koutBC * BC_SPLEEN;
// double death_BC_node   = koutBC * BC_NODE;
// double death_BC_lymph  = koutBC * BC_LYMPH;
// 
// // B-Cells Trafficking Fluxes (Movement)
// double traf_BC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * BC_BLOOD;
// double traf_BC_spleen_node  = ktn * BC_SPLEEN;
// double traf_BC_node_lymph   = knl * BC_NODE;
// double traf_BC_lymph_blood  = klp * BC_LYMPH;

// T-Cells: Production & Death
// We also clamp AF_TC to prevent errors if homeostasis dips slightly below 0
double prod_TC         = kinTC * pow(CLAMP(AF_TC), r);
double death_TC_blood  = koutTC * TC_BLOOD;
double death_TC_spleen = koutTC * TC_SPLEEN;
double death_TC_node   = koutTC * TC_NODE;
double death_TC_lymph  = koutTC * TC_LYMPH;

// T-Cells: Trafficking Fluxes
// Logic: We can only move cells that actually exist
double traf_TC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * TC_BLOOD;
double traf_TC_spleen_node  = ktn * TC_SPLEEN;
double traf_TC_node_lymph   = knl * TC_NODE;
double traf_TC_lymph_blood  = klp * TC_LYMPH;

// B-Cells: Production & Death
double prod_BC         = kinBC * pow(CLAMP(AF_BC), r);
double death_BC_blood  = koutBC * BC_BLOOD;
double death_BC_spleen = koutBC * BC_SPLEEN;
double death_BC_node   = koutBC * BC_NODE;
double death_BC_lymph  = koutBC * BC_LYMPH;

// B-Cells: Trafficking Fluxes
double traf_BC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * BC_BLOOD;
double traf_BC_spleen_node  = ktn * BC_SPLEEN;
double traf_BC_node_lymph   = knl * BC_NODE;
double traf_BC_lymph_blood  = klp * BC_LYMPH;