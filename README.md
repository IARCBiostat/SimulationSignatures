# Simulation Study: Feature Selection and Collider Bias in Causal Interpretation

This repository contains the code for the simulation study presented in **On the construction of molecular signatures of lifestyle exposures**, evaluating feature selection strategies under different latent variable scenarios.

## 📁 Repository Structure

```
simulation_signatures/
├── run_simulation.R                      # Main simulations (Scenarios 1-3)
├── R/
│   ├── analysis_fxns.R                   # Analysis functions
│   ├── generate_data.R                   # Data generation functions
│   ├── generate_figures.R                # Main figure generation
│   ├── generate_figures_supplementary.R  # Supplementary figure generation
│   ├── run_scenario_1.R                  # Script to run scenario 1
│   ├── run_scenario_2.R                  # Script to run scenario 2
│   ├── run_scenario_3.R                  # Script to run scenarios 3
│   ├── run_scenario_2.sh                 # Bash script to run scenario 2
│   ├── run_scenario_3.sh                 # Bash script to run scenario 3
├── config/
│   ├── scenario_1.R                      # Config file for Scenario 1 
│   ├── scenario_2.R                      # Config file for Scenario 2 
│   ├── scenario_3.R                      # Config file for Scenario 3 
├── renv/                                 # Reproducible R environment
├── renv.lock                             # Package versions used in the study
└── README.md                             # You're here
```

## Clone the repository

If you have Git installed, run:
```bash
git clone https://github.com/IARCBiostat/SimulationSignatures/
cd SimulationSignatures
```

Alternatively, you can download the ZIP file directly from GitHub and unzip it.

## Set up the R environment

To ensure the correct package versions, we recommend using the [`renv`](https://rstudio.github.io/renv/) environment:

```r
install.packages("renv")        # If not already installed
renv::restore()                 # Restores package versions from renv.lock
```

This will install the specific versions of all required packages as used in the original analysis.

## Running Simulations

Simulations are designed to be run **in parallel**. The number of cores used is controlled via the `NCore` parameter defined in the corresponding config file (e.g., `config/scenario_1.R`).

### Run main simulations (Scenarios 1-3)

```bash
Rscript run_simulation.R 
```

#### Scenario 2 and 3 (SLURM-based parallelization)

Scenario 2 and 3 is parallelized using SLURM job arrays rather than R's internal parallelization.

Each job processes one parameter combination.

Make sure the array size matches the number of parameter combinations:

```r
nrow(AllParams)
```

Each job uses `SLURM_ARRAY_TASK_ID` to select its task.

### Run supplementary simulations

```bash
Rscript run_supplementary.R config/supplementary_scenario_4.R
```

You can edit the `NCore` value in the config file to match the number of CPU cores you wish to allocate.

### Output

Simulation results will be saved in the `results/` directory as `.rds` files.

## Generating Figures

Once the simulations are complete, you can generate the figures used in the paper and supplementary material:

### Main paper figures

```r
source("R/generate_figures.R")
```

### Supplementary figures

```r
source("R/generate_figures_supplementary.R")
```

Figures will be saved to a `results/Figures/` or equivalent output directory specified in the figure generation script.

## Questions?

For issues, bugs, or questions, feel free to open an issue or contact us.
