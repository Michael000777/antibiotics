# Antibiotics
Scripts to process data from antibiotics experiments.

Used in manuscript https://www.biorxiv.org/content/10.1101/2022.08.04.502802v1


## Contents
- **ExperimentalData**

  All experimental data from plates for the MSI curve measurements (ChangeRHO, ChangeEPS), and for kill curve measurements (KillCurves) in Fig1.
- **MeasureKillCurves**

  Code to generate Fig1.
- **MeasureMSP**

  Wrapper for experimental data (platereaderclass.py) and rendering a figure for plates (plateimageclass.py), together with main fitting code for MSI curves (plates_EstimateMSP.py)
- **Plots**

  Notebook to generate all figures in manuscript.
- **SimulationDynamics**

  Code to numerically solve model equations (PopulationDynamicsClasses.py)
