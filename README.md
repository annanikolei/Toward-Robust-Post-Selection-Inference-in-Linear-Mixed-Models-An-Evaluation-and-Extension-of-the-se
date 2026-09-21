# Repository for Toward Robust Post-Selection Inference in Linear Mixed Models: An Evaluation and Extension of the selfmade Approach

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

## Reproducibility

* The `main_...` scripts in `simulations/` reproduce the simulation studies.
* The `Selfmade/` folder provides a standalone implementation of the proposed method.
* All required result files for evaluation and plotting are included where needed.

## Requirements

* R
* (additional packages may be required depending on the script)

## Notes

* Simulation scripts are structured such that results are generated first, followed by evaluation and plotting.
* Existing result files can be reused to reproduce figures without rerunning simulations.
