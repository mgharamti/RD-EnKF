# RD-EnKF manuscript: El Gharamti, MWR (2022)

This directory contains matlab and shell scripts used 
to generate figures for: 
"A Randomized Dormant Ensemble Kalman Filter" 
by Mohamad El Gharamti
which was submitted to Monthly Weather Review.

## Main Scripts
1. **Scalar 1D Cases:**
    - *gen_Fig_1.m:* A matlab script that generates the Bayesian solution for the RD-EnKF
    - *gen_Figs_3_4.m:* This matlab script runs the 1D scalar example for both linear and nonlinear cases   
2. **Lorenz-96 Attractor:**
    - *run_l96_exps.sh:* A bash shell script that runs all L96 experiments (designed for HPC architecture)
    - *gen_Figs_5_6_7.m:* A matlab script to generate L96 results with perfect modeling conditions
      - Time-series plot with and without adaptive inflation
      - Localization sensitivity 
      - Ensemble correlations plot
    - *gen_Fig_8.m:* Produces L96 results with model errors with and without adaptive inflation
    - *gen_Fig_9.m:* Rank histograms for 3 different experiments
3. **Flood Prediction Using HydroDART:**
    - *gen_Fig_10.m:* Generates the stream network on the flooding domain including gauges and major cities and rivers 
    - *gen_Fig_11.m:* Generates hydrographs comparing EnKF and RD-EnKF at three different gauge locations
    - *gen_Fig_12.m:* Generates CRPSS boxplot summary for all gauges

## Utility Functions:
    - *HydroDARTdiags.m:* An all purpose diagnostic function to HydroDART (generates hydrograph, maps, etc)
    - *gauges_2_indices_subset.m:* A utility function to read gauge IDs and their associated DART indicies
    - *obs_increment_eakf.m:* Computes observation increments for the analysis
    - *move_forward.m:* Integrates the scalar model forward in time
    - *rh.m:* Finds ensemble bins to plot rank histograms

