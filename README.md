# Antibiotics
Scripts to process data from antibiotics experiments.

Manuscript:

Lukas Geyrhofer, Philip Ruelens, Andrew D. Farr, Diego Pesce, J. Arjan G. M. de Visser, Naama Brenner
Minimal Surviving Inoculum in Collective Antibiotic Resistance
mBio, Volume 14  Issue 2  e02456-22, 2023
https://doi.org/10.1128/mbio.02456-22


## Contents
- **ExperimentalData**

  All experimental data from plates for the MSI curve measurements (ChangeRHO, ChangeEPS), and for kill curve measurements (KillCurves) in Fig1.
- **MeasureKillCurves**

  Additional Code to analyze KillCurves, used in Fig1.
- **MeasureMSP**

  Wrapper for experimental data (platereaderclass.py) and rendering a figure for plates (plateimageclass.py), together with main fitting code for MSI curves (plates_EstimateMSP.py)
- **Plots**

  Notebook to generate all figures in manuscript.
- **SimulationDynamics**

  Code to numerically solve model equations (PopulationDynamicsClasses.py)
