// Calculate concentration from amount
// double C_plasma = (PLASMA > 0) ? (PLASMA / Vplasma) : 0.0;
// double C_leaky  = (LEAKY  > 0) ? (LEAKY  / Vleaky)  : 0.0;
// double C_tight  = (TIGHT  > 0) ? (TIGHT  / Vtight)  : 0.0;
// double C_spleen = (SPLEEN > 0) ? (SPLEEN / Vspleen) : 0.0;
// double C_node   = (NODE   > 0) ? (NODE   / Vnode)   : 0.0;
// double C_lymph  = (LYMPH  > 0) ? (LYMPH  / Vlymph)  : 0.0;

double C_plasma = PLASMA / Vplasma;
double C_leaky  = LEAKY  / Vleaky;
double C_tight  = TIGHT  / Vtight;
double C_spleen = SPLEEN / Vspleen;
double C_node   = NODE   / Vnode;
double C_lymph  = LYMPH  / Vlymph;