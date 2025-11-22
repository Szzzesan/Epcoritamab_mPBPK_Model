// TRAFFICKING EQUATIONS
dxdt_INJ = -kdecay * INJ;
dxdt_AF_TC = kt * (1 - AF_TC);
dxdt_AF_BC = kt * (1 - AF_BC);

// T-Cells
double prod_TC  = kinTC * pow(AF_TC, r);
double death_TC_blood  = koutTC * TC_BLOOD;
double death_TC_spleen = koutTC * TC_SPLEEN;
double death_TC_node   = koutTC * TC_NODE;
double death_TC_lymph  = koutTC * TC_LYMPH;

double traf_TC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * TC_BLOOD;
double traf_TC_spleen_node  = ktn * TC_SPLEEN;
double traf_TC_node_lymph   = knl * TC_NODE;
double traf_TC_lymph_blood  = klp * TC_LYMPH;

dxdt_TC_BLOOD  = prod_TC - death_TC_blood - traf_TC_blood_spleen + traf_TC_lymph_blood;
dxdt_TC_SPLEEN = -death_TC_spleen + traf_TC_blood_spleen - traf_TC_spleen_node;
dxdt_TC_NODE   = -death_TC_node + traf_TC_spleen_node - traf_TC_node_lymph;
dxdt_TC_LYMPH  = -death_TC_lymph + traf_TC_node_lymph - traf_TC_lymph_blood;

// B-Cells
double prod_BC  = kinBC * pow(AF_BC, r);
double death_BC_blood  = koutBC * BC_BLOOD;
double death_BC_spleen = koutBC * BC_SPLEEN;
double death_BC_node   = koutBC * BC_NODE;
double death_BC_lymph  = koutBC * BC_LYMPH;

double traf_BC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * BC_BLOOD;
double traf_BC_spleen_node  = ktn * BC_SPLEEN;
double traf_BC_node_lymph   = knl * BC_NODE;
double traf_BC_lymph_blood  = klp * BC_LYMPH;

dxdt_BC_BLOOD  = prod_BC - death_BC_blood - traf_BC_blood_spleen + traf_BC_lymph_blood;
dxdt_BC_SPLEEN = -death_BC_spleen + traf_BC_blood_spleen - traf_BC_spleen_node;
dxdt_BC_NODE   = -death_BC_node + traf_BC_spleen_node - traf_BC_node_lymph;
dxdt_BC_LYMPH  = -death_BC_lymph + traf_BC_node_lymph - traf_BC_lymph_blood;