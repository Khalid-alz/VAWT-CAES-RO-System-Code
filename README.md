# MATLAB Code for Multi-Objective Optimisation of a VAWT-CAES-RO System

This repository contains the MATLAB source code and associated data files to support the findings of the paper: [To be included upon publication].

The code facilitates the simulation, analysis, and multi-objective optimisation of a renewable energy system combining a Vertical Axis Wind Turbine (VAWT), Compressed Air Energy Storage (CAES), and a Reverse Osmosis (RO) desalination unit.

[![DOI](https://zenodo.org/badge/1032837679.svg)](https://doi.org/10.5281/zenodo.16749275)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## System Requirements

* MATLAB R2022b or later.
* Required MATLAB Toolboxes:
    * Parallel Computing Toolbox
    * Optimisation Toolbox
    * Global Optimisation Toolbox
    * Statistics and Machine Learning Toolbox

## Data Description

The analysis relies on several key data files located in the `data/` directory.

* **`ALL_VAWTs_CAES_RO_Cases.xlsx`**: The primary input file, containing all system design configurations from the VAWT-CAES-RO Simulink model, is used for the main simulation and analysis.

* **`simulink_random_LCOW_CO2_filtered.xlsx`**: A validation dataset, generated using Latin Hypercube Sampling from the VAWT-CAES-RO Simulink model, is used to test the accuracy of the three surrogate models (Qp, LCOW, GWP).

* **`RO_WAVE_Grid_Samples.xlsx`**: A dataset containing 125 data points generated using WAVE software, which was used to train the surrogate models for the reverse osmosis unit.

* **`.mat` files**: These files contain either pre-trained surrogate models (`O_GPR_Qp_10F.mat`, `ONN_CO2_10F.mat`, `ONN_LCOW_10F.mat`, `RO_Cp.mat`, etc.) or intermediate data products used by the analysis scripts (`All_data_MOO.mat`).

## Installation and Setup

1.  Clone or download this repository to your local machine.
2.  Open MATLAB.
3.  Set the MATLAB working directory to the root folder of this repository.
4.  Ensure all data files listed above are present in their respective locations.

## Usage: Code Workflow

The scripts are organised to follow the research workflow, from initial system simulation through to validation.

### Phase 1: System Simulation

These scripts run the Simulink model for the VAWT-CAES-RO system.

* **`RUN_Simulink_SINGLE.m`**: Runs the Simulink model for a single, user-defined case. Useful for testing and debugging.
* **`RUN_Simulink_POOL.m`**: Runs the full set of simulations for all cases defined in the input grid file. This script utilises parallel computing to accelerate the process.
* **`RUN_Simulink_RandomSet.m`**: Runs simulations for the random validation set, generated via Latin Hypercube Sampling.

* *Inputs:`Sim_main.slx` The Simulink model for the VAWTs-CAES-RO system is called from the three codes above.


### Phase 2: Economic and Environmental Analysis

This script processes the simulation output to calculate key performance indicators.

* **`LCOW_calculation.m`**: Calculates the Levelised Cost of Water (LCOW) for each case specified in `ALL_VAWTs_CAES_RO_Cases.xlsx`. The results are saved to a new file, `ALL_VAWTs_CAES_RO_Cases_with_LCOW.xlsx`.

### Phase 3: Surrogate Model Validation

These scripts validate the accuracy of the surrogate models against real or test data.

* **`Compare_real_VS_models_Qp.m`**: Compares the permeate flow rate (Qp) surrogate model predictions against test data, calculating the Root Mean Square Error (RMSE) and generating comparison plots.
* **`Compare_real_VS_model_LCOW.m`**: Compares the LCOW surrogate model predictions against test data, calculating RMSE and generating plots.
* **`Compare_real_VS_models_CO2.m`**: Compares the GWP (CO2) surrogate model predictions against test data, calculating RMSE and generating plots.

### Phase 4: Multi-Objective Optimisation

This script performs the core optimisation task using pre-trained surrogate models.

* **`MOO_Qp_LCOW_CO2.m`**: Applies multi-objective optimisation to find the Pareto optimal front for three conflicting objectives: permeate flow rate (Qp), LCOW, and Global Warming Potential (GWP). This script generates the data for **Figure 13**.
    * *Inputs:* `O_GPR_Qp_10F.mat`, `ONN_CO2_10F.mat`, `ONN_LCOW_10F.mat` (surrogate models).

### Phase 5: Multi-Criteria Decision Making

This script analyses the Pareto front solutions to identify a single preferred solution.

* **`TOPSIS.m`**: Applies the Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) to the Pareto optimal solutions. It analyses four scenarios with different objective weightings (as described in **Table 3**) and generates the Parallel Coordinate Plot shown in **Figure 14**.
    * *Input:* `All_data_MOO.mat` file contains data generated from the Multi-Objective Optimization (save its output).


## Licence

This project is licensed under the MIT Licence. See the `LICENCE` file for details.

## How to Cite

If you use this code in your research, please cite both the repository and the associated paper.

1.  **Paper Citation:**
    > [the full paper citation here once it is published]

2.  **Software Citation:**
    > Alzahrani, K. (2025). Code for: A Multi-Objective Optimisation Framework for a VAWT-CAES-RO System (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.16749276