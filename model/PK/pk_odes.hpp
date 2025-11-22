// calculate fluxes & derivatives
// This file assumes concentrations, such as C_plasma are already defined
double dist_plasma_leaky  = Lleaky  * (1.0 - sigma_leaky)  * C_plasma;
double dist_plasma_tight  = Ltight  * (1.0 - sigma_tight)  * C_plasma;
double dist_plasma_spleen = Lspleen * (1.0 - sigma_spleen) * C_plasma;
double dist_leaky_lymph   = Lleaky  * (1.0 - sigma_lymph)  * C_leaky;
double dist_tight_lymph   = Ltight  * (1.0 - sigma_lymph)  * C_tight;
double dist_spleen_node   = Lspleen * (1.0 - sigma_lymph)  * C_spleen;
double dist_node_lymph    = Lspleen * C_node;
double dist_lymph_plasma  = L * C_lymph;

double elim_linear    = CL * C_plasma;
double elim_nonlinear = (Vmax * C_plasma) / (Km + C_plasma) * Vplasma;

dxdt_DEPOT  = -ka * DEPOT;
dxdt_PLASMA = -elim_linear - elim_nonlinear - dist_plasma_leaky - dist_plasma_tight - dist_plasma_spleen + dist_lymph_plasma;
dxdt_LEAKY  = dist_plasma_leaky - dist_leaky_lymph;
dxdt_TIGHT  = dist_plasma_tight - dist_tight_lymph;
dxdt_SPLEEN = dist_plasma_spleen - dist_spleen_node;
dxdt_NODE   = dist_spleen_node - dist_node_lymph;
dxdt_LYMPH  = (ka * DEPOT) + dist_leaky_lymph + dist_tight_lymph + dist_node_lymph - dist_lymph_plasma;