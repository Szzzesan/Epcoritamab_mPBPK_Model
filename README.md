# Epcoritamab mPBPK-QSP Model

## Overview
This repository contains a quantitative systems pharmacology (QSP) model for **Epcoritamab** (CD3xCD20 bispecific antibody), implemented in R and C++ using `mrgsolve`. 

The model implements the work of *Li et al. (2022)*, simulating the interactions between the drug, T-cells, B-cells, and tumor cells across multiple physiological compartments (Blood, Spleen, Lymph Node, Lymph).

## Key Features
* **Modular Design:** Built from 5 integrated submodels (PK, Trafficking, Binding, Activation, Killing).
* **Mechanistic Detail:** Explicitly models Trimer formation ($Drug-CD3-CD20$) as the driver of efficacy.

## Repository Structure
├── model/ 
│ ├── epc_complete.cpp # The MASTER model file  
│ ├── PK/ # (Development) Isolated PK submodel 
│ └── PD/ # (Development) Isolated PD submodels 

├── scripts/ 
│ ├── run_complete.R # MAIN SCRIPT: Runs the full simulation & plots (Run this one) 
│ └── run_tests.R # (Legacy) Unit tests for individual submodels 

├── data/ # not used here

├── output/
│ ├── Epcoritamab_mPBPK_QSP.Rmd # R Markdown source for the final report
│ └── Epcoritamab_mPBPK_QSP.html # html final report

└── Epcoritamab_mPBPK_Model.Rproj

## Getting Started

### Prerequisites
* R (>= 4.0)
* `mrgsolve`
* `tidyverse`

### Running the Model
1.  Clone this repository.
2.  Open `Epcoritamab_mPBPK_Model.Rproj` in RStudio.
3.  Open **`scripts/run_complete.R`**.
4.  Run the script to compile the model and generate plots.

## Model Details
### Mechanism of Action (MoA)
Epcoritamab is a bispecific antibody that simultaneously binds to **CD3** on T-cells and **CD20** on B-cells (or tumor cells).
1.  **Bridge Formation:** The drug physically links a cytotoxic T-cell to a target B-cell, forming a **Trimer** ($Tcell-Drug-Bcell$).
2.  **Activation:** Trimer formation triggers T-cell activation and expansion.
3.  **Killing:** Activated T-cells release cytotoxic granules, lysing the target B-cell (or tumor cell).

---

### Submodel Breakdown

The model is composed of six integrated submodels, defined in `model/epc_complete.cpp`.

#### 1. PK Submodel (Pharmacokinetics)
* **Description:** A two-pore minimal PBPK model describing the distribution of the antibody across physiological compartments (Plasma, Tight Tissue, Leaky Tissue, Spleen, Lymph Node, Lymph).
* **Key Features:**
    * Uses **reflection coefficients ($\sigma$)** to govern convective transport between plasma and lymph.
    * Includes both **Linear Clearance ($CL$)** and **Non-linear Clearance ($V_{max}, K_m$)** (Target-Mediated Drug Disposition).
* **Attention Point:** The model assumes instantaneous mixing within each compartment volume. Distribution is driven by lymph flow rates ($L$).

#### 2. Lymphocyte Trafficking Submodel
* **Description:** Describes the recirculation of T-cells and B-cells between blood, spleen, lymph nodes, and lymph.
* **Key Features:**
    * **Injection Effect (`INJ`):** Simulates the transient margination (drop) of lymphocytes from the blood into the spleen immediately following drug administration.
    * **Homeostasis (`AF_TC`, `AF_BC`):** Adaptive feedback loops sensing total cell counts. If cells are depleted, the system upregulates production (`kin`) to return to baseline.
* **Attention Point:** Correct initialization is critical. The initial cell counts in each tissue must be balanced with trafficking rates ($k_{pt}, k_{tn}$, etc.) to ensure a steady state before dosing.

#### 3. Binding Submodel
* **Description:** Explicit mass-action kinetics for the formation of binary complexes (Dimers) and ternary complexes (Trimers).
* **Key Features:**
    * Tracks **Free**, **Dimer** ($Drug-CD3$, $Drug-CD20$), and **Trimer** species in all lymphoid tissues.
    * **Target Competition:** In the Lymph Node, the drug must choose between binding healthy B-cells or Tumor cells.
    * **Conservation of Mass:** Includes "Recycling Logic"—when a target cell dies (naturally or by killing), the Trimer breaks, and the CD3-bound antibody is released back to the T-cell pool.
* **Attention Point:** This submodel creates a "stiff" ODE system because binding reactions (microseconds) occur much faster than cell growth (days). The parameter `scale_binding` is available to slow down kinetics for numerical stability if needed, though `1.0` (real physics) is preferred.

#### 4. T-Cell Activation Submodel
* **Description:** Links Trimer formation to T-cell proliferation.
* **Key Features:**
    * **Activation Delay (`TAD`):** A biological delay (default 3 days) prevents immediate activation upon dosing.
    * **Signal Transduction:** The rate of activation is proportional to the number of Trimers per target cell.
    * **Expansion:** Activated T-cells (`vATC`) undergo clonal expansion (`pATC`), creating a "cytotoxic army."
* **Attention Point:** The activation parameter `sim_slope` is highly potent. A value too high can cause a "flash fire" effect where targets are wiped out instantly.

#### 5. B-Cell Killing & 6. Tumor Killing Submodels
* **Description:** Describes the elimination of healthy B-cells (pharmacodynamics) and Tumor cells (efficacy).
* **Key Features:**
    * **Fractional Killing:** The killing rate depends on the density of activated T-cells and target cells.
    * **Spherical Geometry:** For the tumor, only cells in the "Reachable Shell" (defined by `depth`) can bind drug and be killed. The necrotic core is protected.
* **Attention Point:**
    * **Units:** The parameters `kkill` and `sim_slope` from literature are often scaled for "Millions of Cells." In this model, they are re-scaled to "Single Cells" (multiplied/divided by $10^6$) to allow for precise tracking of low cell counts.
    * **Consumption Trap:** High killing rates can deplete the "anchor" (B-cells) so fast that Trimer concentrations crash. This is a correct mechanistic result but requires careful parameter tuning.

## References
Li, T. et al. *Semimechanistic Physiologically-Based Pharmacokinetic/Pharmacodynamic Model Informing Epcoritamab Dose Selection...* Clinical Pharmacology & Therapeutics (2022).