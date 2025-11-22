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

$CMT 
// PK Compartments
DEPOT PLASMA LEAKY TIGHT SPLEEN NODE LYMPH
// Trafficking Compartments
TC_BLOOD TC_SPLEEN TC_NODE TC_LYMPH
BC_BLOOD BC_SPLEEN BC_NODE BC_LYMPH
INJ AF_TC AF_BC

$GLOBAL
#define CLAMP(x) ((x) > 0 ? (x) : 0.0)

$ODE @!audit
#include "../../PK/pk_vars.hpp"
#include "../../PK/pk_odes.hpp"
#include "traff_vars.hpp"
#include "traff_odes.hpp"

$TABLE
#include "traff_vars.hpp"

$CAPTURE traf_TC_blood_spleen 
// capture traf_TC_blood_spleen so we can plot 
// to see the trafficking rate spiked 
// when injecting the drug