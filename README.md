# MATLAB code for ZOI model figures

This repository contains the MATLAB scripts used to generate the computational figures for the accompanying paper.

Each script is self-contained and can be run independently. Generated figures are saved to the `outputs` directory.

## Files

| File | Figure |
|---|---:|
| `figure_2_1.m` | Figures 2.1a–b |
| `figure_2_2.m` | Figure 2.2 |
| `figure_2_3.m` | Figure 2.3 |
| `figure_3_1.m` | Figures 3.1a–d |
| `figure_3_2.m` | Figures 3.2a–b |
| `figure_3_3.m` | Figures 3.3a–b |
| `figure_4_1.m` | Figures 4.1a-b |
| `figure_4_2.m` | Figures 4.2a–b |
| `figure_D_1.m` | Figures D.1a–b |

## Requirements

The scripts require MATLAB with support for local functions in scripts.

The following toolboxes are required:

- **Optimization Toolbox** for `fsolve`;
- **Communications Toolbox** for `marcumq`;
- **Statistics and Machine Learning Toolbox** for functions including `prctile` and random sampling used in the sensitivity analysis.

## Running the code

Open MATLAB, set the current folder to the repository directory, and run the desired script. For example:

```matlab
figure_3_1
```

or

```matlab
run('figure_3_1.m')
```

The required `outputs` directory is created automatically if it does not already exist.

Scripts containing multiple manuscript panels export each panel as a separate PDF file.

## Citation

Please cite the accompanying paper when using this code:

> Patterson, A. C., Bradshaw-Hajek, B. H., Murphy, R. J., Bunder, J. E., Venn, X. L., and Tam, A. K. Y.  
> *Revisiting the zone of inhibition experiments for antimicrobial susceptibility: a diffusion-limited hole opening problem*.

## Contact

Alexander K. Y. Tam  
alexander.tam@adelaide.edu.au
