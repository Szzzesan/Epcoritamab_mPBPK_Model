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

$CMT DEPOT PLASMA LEAKY TIGHT SPLEEN NODE LYMPH

$GLOBAL
#define CLAMP(x) ((x) > 0 ? (x) : 0.0)

$ODE @!audit
// The Engine needs BOTH variables and equations
#include "pk_vars.hpp"
#include "pk_odes.hpp"

$TABLE
// The Reporter needs just the variables to capture C_plasma
#include "pk_vars.hpp"

$CAPTURE C_plasma