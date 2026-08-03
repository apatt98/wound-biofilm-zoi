# wound-biofilm-zoi
Contains the MATLAB scripts used to generate the figures in:

**Revisiting the zone of inhibition experiments for antimicrobial susceptibility: a diffusion-limited hole opening problem**


Each script is self-contained and corresponds to a numbered figure in the paper.

## Abstract

...

## Files

| File | Paper Figure | Description |
|:---:|:---:|:---|
| 'figure_2_1.m' | Figure 2.1 | Antimicrobial concentration profiles for the infinite-domain Dirac initial condition at selected times. |
| 'figure_2_2.m' | Figure 2.2 | Comparison of the MIC level curves for the four combinations of initial condition and domain size. |
| 'figure_2_3.m' | Figure 2.3 | Finite-domain level curves for several domain radii, for both Dirac and Heaviside initial conditions. |
| 'figure_3_1.m' | Figure 3.1 | Bacterial density profiles at selected times for fixed bacterial death parameter. |
| 'figure_3_2.m' | Figure 3.1 | Bacterial density profiles at for several values of the bacterial death parameter. |
| 'figure_4_1.m' | Figure 4.1 | One at a time sensitivity analysis for the meropenem data. |
| 'figure_4_2.m' | Figure 4.2 | Model fit and uncertainty band for the meropenem data. |
| 'figure_4_3.m' | Figure 4.3 | One at a time sensitivity analysis for the ciprofloxacin data. |
| 'figure_4_4.m' | Figure 4.4 | Model fit and uncertainty band for the ciprofloxacin data. |

The 'outputs' folder is used for the generated PDF files.

## Requirements

The scripts require MATLAB R2016b or later. The following toolboxes are required:

- **Communications Toolbox** for 'marcumq' in Figures 3.1 and 3.2
- **Optimization Toolbox** for 'fsolve' in Figures 3.1 and 3.2
- **Statistics and Machine Learning Toolbox** for 'unifrns' and 'prctile' in Figures 4.1--4.4

## Data and modelling notes

The experimental values used in Figures 4.1--4.4 are included directly in the relevant scripts. They were transcribed from the studies cited in the manuscript.

## Citation

...
