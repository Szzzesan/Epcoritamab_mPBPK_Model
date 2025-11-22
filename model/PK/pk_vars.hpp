// Calculate concentration from amount
// double C_plasma = PLASMA / Vplasma;
// double C_leaky  = LEAKY  / Vleaky;
// double C_tight  = TIGHT  / Vtight;
// double C_spleen = SPLEEN / Vspleen;
// double C_node   = NODE   / Vnode;
// double C_lymph  = LYMPH  / Vlymph;

// double C_plasma = (PLASMA > 0) ? (PLASMA / Vplasma) : 0.0;
// double C_leaky  = (LEAKY  > 0) ? (LEAKY  / Vleaky)  : 0.0;
// double C_tight  = (TIGHT  > 0) ? (TIGHT  / Vtight)  : 0.0;
// double C_spleen = (SPLEEN > 0) ? (SPLEEN / Vspleen) : 0.0;
// double C_node   = (NODE   > 0) ? (NODE   / Vnode)   : 0.0;
// double C_lymph  = (LYMPH  > 0) ? (LYMPH  / Vlymph)  : 0.0;

double C_plasma = CLAMP(PLASMA) / Vplasma;
double C_leaky  = CLAMP(LEAKY)  / Vleaky;
double C_tight  = CLAMP(TIGHT)  / Vtight;
double C_spleen = CLAMP(SPLEEN) / Vspleen;
double C_node   = CLAMP(NODE)   / Vnode;
double C_lymph  = CLAMP(LYMPH)  / Vlymph;