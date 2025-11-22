$PARAM @annotated
// --- PK Submodel Parameters (from Table S1) ---
ka          : 0.131   : Absorption rate constant (1/day)               
CL          : 2.47    : Linear clearance (L/day)                    
Km          : 0.0461  : MM constant (nM)                          
Vmax        : 0.185   : Max nonlinear clearance rate (nM/day)   
   
// Volumes (L) 
Vplasma     : 2.6     : Plasma volume                              
Vleaky      : 4.37    : Leaky tissue volume                    
Vtight      : 8.11    : Tight tissue volume                       
Vspleen     : 0.0433  : Spleen volume                         
Vnode       : 0.0082  : Lymph node volume                      
Vlymph      : 5.2     : Lymph volume                               
 
// Lymph Flows (L/day) 
Lleaky      : 1.94    : Flow through leaky tissue              
Ltight      : 0.957   : Flow through tight tissue                
Lspleen     : 0.304   : Flow through spleen                           
L           : 2.9     : Total lymph flow                          
 
// Reflection Coefficients 
sigma_leaky : 0.85    : Reflection coeff for leaky tissue                          
sigma_tight : 0.95    : Reflection coeff for tight tissue                    
sigma_spleen: 0.85    : Reflection coeff for spleen                
sigma_lymph : 0.2     : Reflection coeff for lymph                         

// --- Lymphocyte Trafficking Parameters (from Table S1) ---
// T-Cell Parameters
kinTC       : 310e6   : Production rate of T cells (cells/day) 
koutTC      : 0.00422 : Natural death rate of T cells (1/day) 
TC_base     : 1469    : Baseline T-cell count in blood (cells/mm3) 

// B-Cell Parameters
kinBC       : 223e6   : Production rate of B cells (cells/day) 
koutBC      : 0.0189  : Natural death rate of B cells (1/day) 
BC_base     : 236     : Baseline B-cell count in blood (cells/mm3) 

// Trafficking Rates (1/day)
kpt         : 50      : Blood to Spleen 
ktn         : 1.67    : Spleen to Node 
knl         : 1672    : Node to Lymph 
klp         : 2.63    : Lymph to Blood 

// Injection Effect (Transient drop in lymphocytes)
kdecay      : 0.681   : Duration of injection effect (1/day) 
INJ_Scaler  : 9.08    : Magnitude of injection effect 

// Homeostasis
kt          : 0.00027 : Lymphocyte homeostasis control rate (1/day) 
r           : 2.24    : Strength of lymphocyte homeostasis control 

// --- BINDING PARAMETERS --- (From Table S1)
konCD3      : 18.1    : Assoc rate to CD3 (L/nmol/day)
koffCD3     : 285     : Dissoc rate from CD3 (1/day)
konCD20     : 4.15    : Assoc rate to CD20 (L/nmol/day)
koffCD20    : 22.5    : Dissoc rate from CD20 (1/day)
kdegCD3     : 1.584   : Degradation of free CD3 (1/day)
kdegCD20    : 1.584   : Degradation of free CD20 (1/day)
kintCD3     : 1.584   : Internalization of Ab-CD3 dimer (1/day)
kintCD20    : 1.584   : Internalization of Ab-CD20 dimer (1/day)
R_CD3       : 30000   : CD3 receptors per T-cell
R_CD20      : 100000  : CD20 receptors per B-cell

$CMT
// PK Compartments
DEPOT PLASMA LEAKY TIGHT SPLEEN NODE LYMPH

// Lymphocyte Compartments (TC = T-cells, BC = B-cells)
TC_BLOOD TC_SPLEEN TC_NODE TC_LYMPH
BC_BLOOD BC_SPLEEN BC_NODE BC_LYMPH

// Helper Compartments
INJ     // Injection effect tracker
AF_TC   // Adaptive Feedback for T-cells
AF_BC   // Adaptive Feedback for B-cells

// Binding Compartments (All tissues)
// BLOOD
FREE_CD3_BLOOD  FREE_CD20_BLOOD  DIMER_CD3_BLOOD  DIMER_CD20_BLOOD  TRIMER_BLOOD
// SPLEEN
FREE_CD3_SPLEEN FREE_CD20_SPLEEN DIMER_CD3_SPLEEN DIMER_CD20_SPLEEN TRIMER_SPLEEN
// NODE
FREE_CD3_NODE   FREE_CD20_NODE   DIMER_CD3_NODE   DIMER_CD20_NODE   TRIMER_NODE
// LYMPH
FREE_CD3_LYMPH  FREE_CD20_LYMPH  DIMER_CD3_LYMPH  DIMER_CD20_LYMPH  TRIMER_LYMPH

$MAIN
// Constants
double AVOG = 6.022e23;
double nmol_per_molecule = 1e9 / AVOG; // Helper to convert counts to nmol

// Initial Conditions for Lymphocytes (Optional: Set baselines here if needed)
// For now, we assume the system starts at steady state or 0 and builds up.

// --- AUTO-INITIALIZATION for the binding model ---
// If the time is 0, calculate the starting receptor amounts based on cell counts

if(NEWIND <= 1) {
  // Blood
  FREE_CD3_BLOOD_0  = TC_BLOOD_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_BLOOD_0 = BC_BLOOD_0 * R_CD20 * nmol_per_molecule;
  
  // Spleen
  FREE_CD3_SPLEEN_0  = TC_SPLEEN_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_SPLEEN_0 = BC_SPLEEN_0 * R_CD20 * nmol_per_molecule;
  
  // Node
  FREE_CD3_NODE_0  = TC_NODE_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_NODE_0 = BC_NODE_0 * R_CD20 * nmol_per_molecule;
  
  // Lymph
  FREE_CD3_LYMPH_0  = TC_LYMPH_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_LYMPH_0 = BC_LYMPH_0 * R_CD20 * nmol_per_molecule;
}

$ODE
// ------------------------------------------------------------------
// 1. PK SUBMODEL EQUATIONS
// ------------------------------------------------------------------
double C_plasma = PLASMA / Vplasma;
double C_leaky  = LEAKY  / Vleaky;
double C_tight  = TIGHT  / Vtight;
double C_spleen = SPLEEN / Vspleen;
double C_node   = NODE   / Vnode;
double C_lymph  = LYMPH  / Vlymph;

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

// ------------------------------------------------------------------
// 2. LYMPHOCYTE TRAFFICKING SUBMODEL
// ------------------------------------------------------------------
// The "Injection Effect": Drives cells out of blood into spleen temporarily
dxdt_INJ = -kdecay * INJ;

// Adaptive Feedback (Homeostasis) equations
// We need the current TC_BLOOD count vs baseline to drive production
// Note: Simplified version of Eq 1344 for now (assuming steady state baseline is normalized)
dxdt_AF_TC = kt * (1 - AF_TC); // Placeholder: Full equation requires linking to TC count
dxdt_AF_BC = kt * (1 - AF_BC);

// --- T-CELL TRAFFICKING ---
// Production & Death
double prod_TC  = kinTC * pow(AF_TC, r);
double death_TC_blood = koutTC * TC_BLOOD;

//only restore these after the test of binding succeeds
// double death_TC_spleen = koutTC * TC_SPLEEN;
// double death_TC_node = koutTC * TC_NODE;
// double death_TC_lymph = koutTC * TC_LYMPH;

// Movement Flows
// Key Logic: Blood -> Spleen is boosted by INJ effect
double traf_TC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * TC_BLOOD;
//double traf_TC_spleen_node  = ktn * TC_SPLEEN;
//double traf_TC_node_lymph   = knl * TC_NODE;
double traf_TC_lymph_blood  = klp * TC_LYMPH;

// ODEs for T-Cells
dxdt_TC_BLOOD  = prod_TC - death_TC_blood - traf_TC_blood_spleen + traf_TC_lymph_blood;
dxdt_TC_SPLEEN = traf_TC_blood_spleen;
dxdt_TC_NODE = 0;
dxdt_TC_LYMPH = -traf_TC_lymph_blood;

//dxdt_TC_SPLEEN = -death_TC_spleen + traf_TC_blood_spleen - traf_TC_spleen_node;
//dxdt_TC_NODE   = -death_TC_node + traf_TC_spleen_node - traf_TC_node_lymph;
//dxdt_TC_LYMPH  = -death_TC_lymph + traf_TC_node_lymph - traf_TC_lymph_blood;

// --- B-CELL TRAFFICKING ---
// Production & Death
double prod_BC  = kinBC * pow(AF_BC, r);
double death_BC_blood = koutBC * BC_BLOOD;

//only restore these after the test of binding succeeds
// double death_BC_spleen = koutBC * BC_SPLEEN;
// double death_BC_node = koutBC * BC_NODE;
// double death_BC_lymph = koutBC * BC_LYMPH;

// Movement Flows
double traf_BC_blood_spleen = kpt * (1 + INJ_Scaler * INJ) * BC_BLOOD;
//double traf_BC_spleen_node  = ktn * BC_SPLEEN;
//double traf_BC_node_lymph   = knl * BC_NODE;
double traf_BC_lymph_blood  = klp * BC_LYMPH;

// ODEs for B-Cells
dxdt_BC_BLOOD  = prod_BC - death_BC_blood - traf_BC_blood_spleen + traf_BC_lymph_blood;
dxdt_BC_SPLEEN = traf_BC_blood_spleen; 
dxdt_BC_NODE = 0; 
dxdt_BC_LYMPH = -traf_BC_lymph_blood;
//dxdt_BC_SPLEEN = -death_BC_spleen + traf_BC_blood_spleen - traf_BC_spleen_node;
//dxdt_BC_NODE   = -death_BC_node + traf_BC_spleen_node - traf_BC_node_lymph;
//dxdt_BC_LYMPH  = -death_BC_lymph + traf_BC_node_lymph - traf_BC_lymph_blood;

// --- 3. BINDING SUBMODEL (Blood Only) ---

// Calculate Concentrations of Receptors (nmol/L)
// Amount (nmol) = Count (cells) * Receptors/Cell * (1e9/Avogadro)
// Conc (nM)     = Amount / Volume

// Convert compartment amounts (nmol) to concentrations (nM) for the reaction
double C_FreeCD3   = FREE_CD3_BLOOD / Vplasma;
double C_FreeCD20  = FREE_CD20_BLOOD / Vplasma;
double C_DimerCD3  = DIMER_CD3_BLOOD / Vplasma;
double C_DimerCD20 = DIMER_CD20_BLOOD / Vplasma;
double C_Trimer    = TRIMER_BLOOD / Vplasma;

// Reaction Velocities (Law of Mass Action)
// 1. Formation of Dimers
double v_form_DimerCD3  = konCD3  * C_plasma * C_FreeCD3;
double v_form_DimerCD20 = konCD20 * C_plasma * C_FreeCD20;

// 2. Dissociation of Dimers
double v_break_DimerCD3  = koffCD3  * C_DimerCD3;
double v_break_DimerCD20 = koffCD20 * C_DimerCD20;

// 3. Formation of Trimer (Cross-linking)
double v_form_Trimer_viaCD3  = konCD20 * C_DimerCD3  * C_FreeCD20;
double v_form_Trimer_viaCD20 = konCD3  * C_DimerCD20 * C_FreeCD3;

// 4. Dissociation of Trimer
double v_break_Trimer_toCD3  = koffCD20 * C_Trimer; // Breaks the CD20 bond
double v_break_Trimer_toCD20 = koffCD3  * C_Trimer; // Breaks the CD3 bond

// Receptor Synthesis & Degradation
// Synthesis Rate = kdeg * R_per_cell * Cell_Count (converted to nmol)
double syn_CD3  = kdegCD3  * TC_BLOOD * R_CD3  * nmol_per_molecule;
double syn_CD20 = kdegCD20 * BC_BLOOD * R_CD20 * nmol_per_molecule;

// Differential Equations (The Changes)
// Note: We multiply rates by Vplasma to convert Conc/Time -> Amount/Time

dxdt_FREE_CD3_BLOOD = ( -v_form_DimerCD3 + v_break_DimerCD3
                        -v_form_Trimer_viaCD20 + v_break_Trimer_toCD20
                        -kdegCD3 * C_FreeCD3 * Vplasma
                        + syn_CD3 );
                          
dxdt_FREE_CD20_BLOOD = ( -v_form_DimerCD20 + v_break_DimerCD20
                         -v_form_Trimer_viaCD3 + v_break_Trimer_toCD3
                         -kdegCD20 * C_FreeCD20 * Vplasma
                         + syn_CD20 );
                           
dxdt_DIMER_CD3_BLOOD = ( v_form_DimerCD3 - v_break_DimerCD3
                        -v_form_Trimer_viaCD3 + v_break_Trimer_toCD3 // Trimer formation consumes dimer
                        +v_break_Trimer_toCD20 // Trimer break releases dimer
                        -kintCD3 * C_DimerCD3 * Vplasma );
                                                      
dxdt_DIMER_CD20_BLOOD = ( v_form_DimerCD20 - v_break_DimerCD20
                          -v_form_Trimer_viaCD20 + v_break_Trimer_toCD20
                          +v_break_Trimer_toCD3
                          -kintCD20 * C_DimerCD20 * Vplasma );
                                                                                  
dxdt_TRIMER_BLOOD     = ( v_form_Trimer_viaCD3 + v_form_Trimer_viaCD20
                          -v_break_Trimer_toCD3 - v_break_Trimer_toCD20 );

// Update Plasma Drug Loss due to binding (Binding consumes drug)
// We subtract the net drug used to form dimers
dxdt_PLASMA = dxdt_PLASMA - (v_form_DimerCD3 - v_break_DimerCD3)*Vplasma
- (v_form_DimerCD20 - v_break_DimerCD20)*Vplasma;

$CAPTURE
C_plasma C_Trimer