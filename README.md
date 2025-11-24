# Epcoritamab mPBPK-QSP Model

## Overview
This repository contains a quantitative systems pharmacology (QSP) model for **Epcoritamab** (CD3xCD20 bispecific antibody), implemented in R and C++ using `mrgsolve`. 

The model replicates the work of *Li et al. (2022)*, simulating the interactions between the drug, T-cells, B-cells, and tumor cells across multiple physiological compartments (Blood, Spleen, Lymph Node, Lymph).

## Key Features
* **Modular Design:** Built from 5 integrated submodels (PK, Trafficking, Binding, Activation, Killing).
* **Mechanistic Detail:** Explicitly models Trimer formation ($Drug-CD3-CD20$) as the driver of efficacy.
* **Numerical Stability:** Handles stiff ODE systems (fast binding vs. slow growth) using robust scaling and solver settings.

## Repository Structure
├── model/ 
│ ├── epc_complete.cpp # The MASTER model file (Run this one) 
│ ├── PK/ # (Development) Isolated PK submodel 
│ └── PD/ # (Development) Isolated PD submodels 

├── scripts/ 
│ ├── run_complete.R # MAIN SCRIPT: Runs the full simulation & plots 
│ ├── run_tests.R # (Legacy) Unit tests for individual submodels 
│ └── Epcoritamab_Report.Rmd # R Markdown source for the final report 

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
4.  Run the script to compile the model and generate the diagnostic dashboard.

## Model Details
The model solves a system of ~50 differential equations. Key mechanisms include:
* **T-Cell Activation Delay (TAD):** Activation logic initiates 3 days post-dose.
* **Homeostasis:** Adaptive feedback loops maintain lymphocyte counts in the absence of drug.
* **Fractional Killing:** Cell death rates are coupled to receptor/trimer loss to ensure mass conservation.

## References
Li, T. et al. *Semimechanistic Physiologically-Based Pharmacokinetic/Pharmacodynamic Model Informing Epcoritamab Dose Selection...* Clinical Pharmacology & Therapeutics (2022).