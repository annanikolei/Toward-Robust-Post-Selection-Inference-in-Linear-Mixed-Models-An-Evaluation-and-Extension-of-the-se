# Repository for: Toward Robust Post-Selection Inference in Linear Mixed Models: An Evaluation and Extension of the selfmade Approach

This repository contains the code accompanying our paper on **valid post-selection inference in linear mixed models**.

## Structure

* `Selfmade/`
  Implementation of the proposed *selfmade* approach (Rügamer et al., 2022).
  Contains a standalone script to run the method.

* `simulations/`
  Code and results for **Simulation Studies 1–4**, including additional and follow-up simulations.
  The `main_...` scripts reproduce the simulations and generate the corresponding plots.

* `illustrative_figures/` *(optional, if you use it)*
  Code for generating conceptual/illustrative figures used in the paper.

* Supplementary Material. The supplementary material is provided in the file Supplement.pdf

## Reproducibility

* The `main_...` scripts in `simulations/` reproduce the simulation studies.
* The `Selfmade/` folder provides a standalone implementation of the proposed method.
* All required result files for evaluation and plotting are included where needed.

## Requirements

* R
* (additional packages may be required depending on the script)

## Before running a script

Set the repository root (the folder containing this README) as R's working directory.

The main scripts are intended to be started from the repository root. For example:

`source("simulations/main_SimStudy1_Results.R")`

## Prepare the supplied result files

The ZIP archives in `simulations/` contain result files. Before running a main results script, extract the relevant archive into a folder with the same name as the archive, without `.zip`. For example, extract `SimStudy1_Results.zip` so that the files are in:

`
simulations/SimStudy1_Results/
`

Similarly, extract `SimulationResults_supplied_D.zip` into `simulations/SimulationResults_supplied_D/`. 

## Reproduce simulations and figures

The main script `simulations/main_SimStudy4_Results_originalSimulation.R` has two modes, controlled by `RUN_SIMULATIONS` near the top of the file. Start R with the repository root as the working directory and source the script, for example:

`
source("simulations/main_SimStudy1_Results.R")
`

- **Create figures from supplied results:** leave `RUN_SIMULATIONS <- FALSE`. First extract `SimStudy1_Results.zip` into `simulations/SimStudy1_Results/` and `SimulationResults_supplied_D.zip` into `simulations/SimulationResults_supplied_D/`.
- **Run the simulations and then create figures:** set `RUN_SIMULATIONS <- TRUE`. The script creates `simulations/SimStudy1_Results/` if needed and saves simulation results there. The supplied-D results must still be extracted as described above.

## Simulation scripts

Each study's supporting scripts divide the work into data generation (`*_dgp.R`), one simulation (`*_one_simulation.R`), and repeated simulation (`*_simulation_study.R`). The `main_...` scripts connect those components and produce figures. 

