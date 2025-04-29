# Effect of Spatial Heterogeneity in Spatial Metapopulation Models

This repository contains the code and data used in the study:  
**_"Effect of Spatial Heterogeneity in Spatial Metapopulation Models."_**

We explore how small fractions of uninhabitable patches in a landscape can drive qualitative shifts in metapopulation dynamics.  
All simulations were implemented in MATLAB and primarily run on UC Davis high-performance computing clusters (FARM).

---

## 📂 Repository Overview

### 1. `all_analysis_misc`
- **Generates** `all_info.mat`, containing 1920 entries. Each entry includes:
  - Parasitoid start location and density
  - Indices for 1%, 5%, and 10% uninhabitable patches (in a 128×128 lattice)
  - Habitat site classifications (2, 3, or 4 habitat neighbors) across different levels of heterogeneity (0%, 1%, 5%, 10%)
- **Selects** 100 random habitat pairs for temporal correlation calculations.
- **Creates** movies showing the time-lapse evolution of host/parasitoid heatmaps.

### 2. `Perc_death_ab`, `Perc_ref_ref`, `Perc_periodic`
- **Simulate** time series for host density.
- **Calculate** the mean host density over the last 1000 timesteps.
- **Analyze** sites with 2, 3, and 4 habitat neighbors.

### 3. `Perc_death_ab_tempcorrhost`, `Perc_ref_ref_tempcorrhost`
- **Calculate** temporal correlations of host densities across 100 adjacent habitat sites.
- **Focus**: Final 1000 timesteps, and sites with varying habitat connectivity (2, 3, or 4 neighbors).

### 4. `Perc_death_ab_lattice`, `Perc_ref_ref_lattice`
- **Simulate** the entire lattice evolution over the last 1000 timesteps.
- **Output**: Full spatial timeseries used for generating heatmaps and analyzing spatial pattern formation.

### 5. `Ph_OG_death`
- **Example** batch script used to set up and run simulations efficiently on UC Davis FARM computing cluster.

### 6. `par_save`
- **Utility** function to dynamically and safely save output files after simulation runs.

---

## 🛠 Technical Details
- **Language**: MATLAB
- **Computing**: Simulations were run using parallelization on UC Davis HPC clusters when needed.
- **Data Size**: Simulation outputs can be large due to storing full 128×128 lattices over extended timescales.


## 📬 Contact
For questions or collaborations, feel free to reach out via GitHub or email - akushal@ucdavis.edu.

