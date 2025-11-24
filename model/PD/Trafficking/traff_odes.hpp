// Trafficking ODEs

// Injection & Homeostasis
// We need the current TC_BLOOD count vs baseline to drive production
// Note: Simplified version of Eq 1344 for now (assuming steady state baseline is normalized)
dxdt_INJ = -kdecay * INJ;
dxdt_AF_TC = kt * (1 - AF_TC); // Placeholder: Full equation requires linking to TC count
dxdt_AF_BC = kt * (1 - AF_BC);

//T-Cell Derivatives
dxdt_TC_BLOOD  = prod_TC - death_TC_blood - traf_TC_blood_spleen + traf_TC_lymph_blood;
dxdt_TC_SPLEEN = -death_TC_spleen + traf_TC_blood_spleen - traf_TC_spleen_node;
dxdt_TC_NODE   = -death_TC_node + traf_TC_spleen_node - traf_TC_node_lymph;
dxdt_TC_LYMPH  = -death_TC_lymph + traf_TC_node_lymph - traf_TC_lymph_blood;

//B-Cell Derivatives
dxdt_BC_BLOOD  = prod_BC - death_BC_blood - traf_BC_blood_spleen + traf_BC_lymph_blood;
dxdt_BC_SPLEEN = -death_BC_spleen + traf_BC_blood_spleen - traf_BC_spleen_node;
dxdt_BC_NODE   = -death_BC_node + traf_BC_spleen_node - traf_BC_node_lymph;
dxdt_BC_LYMPH  = -death_BC_lymph + traf_BC_node_lymph - traf_BC_lymph_blood;