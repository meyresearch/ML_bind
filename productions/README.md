# productions

This directory contains folders corresponding to each protein-ligand complex simulation. The folders are named according to the ligand used in the simulation, as all simulations have been run using Tyk2 as the protein. The folders containing 'ANI' in the name are for hybrid ANI-2x/AMBER simulations, whereas those without it are for traditional AMBER simulations.

Each folder contains the following files:
- Raw structure files for Tyk2 (.pdb) and the ligand (.sdf)
- AMBER prep file (.prepi) for the ligand
- Processed structure file for the ligand (.mol2)
- Ligand parameter file (.frcmod)
- Topology (.prmtop) and input coordinate (.inpcrd) files for vacuum Tyk2, vacuum ligand, vacuum complex, and solvated complex
- Text files containing the reported system conditions (e.g. temperature, energies, etc.) throughout the MD simulations, for both the equilibration period (equilibration.csv) and production period (production.csv)
- A folder called 'equilibration_graphs' containing plots of various variables throughout the equilibration period of simulations. Plots are saved as images (.png) according to the dependent variable plotted.
- Structure file of the final equilibration simulation frame (equilibration.pdb) and production simulation frame (production.pdb)
- Result from MM/PBSA binding affinity calculations (MMPBSA_results.dat)
- Ligand RMSD results in text format (rmsd_results.csv), as well as plotted and saved as an image (rmsd.png)
- OpenEye protein-ligand interaction results, saved as a vector file (.svg)
- Percentage of frames in which specific protein-ligand interactions are present (results from interactions_frames.py), saved as a .csv file (interactions_frames.csv), as well as plotted and saved as an image (interactions_frames.png)
