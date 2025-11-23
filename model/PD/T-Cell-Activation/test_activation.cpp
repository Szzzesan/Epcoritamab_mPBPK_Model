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
kinTC       : 310e6   : Production rate of T cells (cells/day)
koutTC      : 0.00422 : Natural death rate of T cells (1/day)
kinBC       : 223e6   : Production rate of B cells (cells/day)
koutBC      : 0.0189  : Natural death rate of B cells (1/day)

// Trafficking Rates (1/day)
kpt         : 50      : Blood to Spleen
ktn         : 1.67    : Spleen to Node
knl         : 1672    : Node to Lymph
klp         : 2.63    : Lymph to Blood

// Injection Effect
kdecay      : 0.681   : Duration of injection effect (1/day)
INJ_Scaler  : 9.08    : Magnitude of injection effect

// Homeostasis
kt          : 0.00027 : Lymphocyte homeostasis control rate (1/day)
r           : 2.24    : Strength of lymphocyte homeostasis control

// --- Binding Parameters --- (From Table S1)
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
scale_binding : 0.01  : Scaling factor for stiffness control

// --- Activation Params ---
sim_slope        : 0.007 : Rate of T-cell activation against B cells
sim_slopetumor   : 7e-6  : Rate of T-cell activation against tumor cells
expand_factor    : 9.27  : Clonal expansion rate (1/day)
koutATC          : 0.05  : Elimination of activated T cells (1/day)
TAD              : 3.0   : T-cell activation delay (days)
Tp               : 10.8  : Proliferation duration (days)
Trimer_Threshold : 24    : Trimers/cell needed for full activation

$CMT
// PK
DEPOT PLASMA LEAKY TIGHT SPLEEN NODE LYMPH
// TRAFFICKING
TC_BLOOD TC_SPLEEN TC_NODE TC_LYMPH
BC_BLOOD BC_SPLEEN BC_NODE BC_LYMPH
INJ AF_TC AF_BC
// BINDING
FREE_CD3_BLOOD  FREE_CD20_BLOOD  DIMER_CD3_BLOOD  DIMER_CD20_BLOOD  TRIMER_BLOOD
FREE_CD3_SPLEEN FREE_CD20_SPLEEN DIMER_CD3_SPLEEN DIMER_CD20_SPLEEN TRIMER_SPLEEN
FREE_CD3_NODE   FREE_CD20_NODE   DIMER_CD3_NODE   DIMER_CD20_NODE   TRIMER_NODE
FREE_CD3_LYMPH  FREE_CD20_LYMPH  DIMER_CD3_LYMPH  DIMER_CD20_LYMPH  TRIMER_LYMPH
// ACTIVATION
vATC_BC_BLOOD vATC_BC_SPLEEN vATC_BC_NODE vATC_BC_LYMPH
vATC_TUMOR_NODE
pATC_BLOOD pATC_SPLEEN pATC_NODE pATC_LYMPH

$GLOBAL
#define CLAMP(x) ((x) > 0 ? (x) : 0.0)
double AVOG = 6.022e23;
double nmol_per_molecule = 1e9 / AVOG;

$MAIN
// Initialize Activation Compartments to 0 (Start inactive)
if(NEWIND <= 1) {
  // Binding (Initialize based on baseline cells)
  FREE_CD3_BLOOD_0  = TC_BLOOD_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_BLOOD_0 = BC_BLOOD_0 * R_CD20 * nmol_per_molecule;
  FREE_CD3_SPLEEN_0  = TC_SPLEEN_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_SPLEEN_0 = BC_SPLEEN_0 * R_CD20 * nmol_per_molecule;
  FREE_CD3_NODE_0  = TC_NODE_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_NODE_0 = BC_NODE_0 * R_CD20 * nmol_per_molecule;
  FREE_CD3_LYMPH_0  = TC_LYMPH_0 * R_CD3  * nmol_per_molecule;
  FREE_CD20_LYMPH_0 = BC_LYMPH_0 * R_CD20 * nmol_per_molecule;
  
  // Activation (Start at 0)
  vATC_BC_BLOOD_0 = 0; 
  vATC_BC_SPLEEN_0 = 0; 
  vATC_BC_NODE_0 = 0; 
  vATC_BC_LYMPH_0 = 0;
  vATC_TUMOR_NODE_0 = 0;
  pATC_BLOOD_0 = 0; 
  pATC_SPLEEN_0 = 0; 
  pATC_NODE_0 = 0; 
  pATC_LYMPH_0 = 0;
}

$ODE @!audit
#include "../../PK/pk_vars.hpp"
#include "../../PK/pk_odes.hpp"
#include "../Trafficking/traff_vars.hpp"
#include "../Trafficking/traff_odes.hpp"
#include "../Binding/binding_vars.hpp"
#include "../Binding/binding_odes.hpp"

double TimeAfterDose = TIME;

#include "activation_vars.hpp"
#include "activation_odes.hpp"

$TABLE
#include "../../PK/pk_vars.hpp"
#include "../Trafficking/traff_vars.hpp"
#include "../Binding/binding_vars.hpp"

double TimeAfterDose = TIME;

#include "activation_vars.hpp"
double Total_pATC = CLAMP(pATC_BLOOD) + CLAMP(pATC_SPLEEN) + CLAMP(pATC_NODE) + CLAMP(pATC_LYMPH);

$CAPTURE Total_pATC