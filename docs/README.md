# docs

This directory contains documents with results of the study that are applicable to all protein-ligand simulations.

If you are looking for results that are specific to a certain ligand (e.g. protein-ligand interactions, specific ligand RMSD plot, etc.), please visit the relevant subfolder for the ligand in the 'productions' directory.

This directory contains the following files:
- MMPBSA_results.csv : contains MM/PBSA binding affinity results for each protein-ligand system, as well as experimental binding affinity results obtained from Wang et al.
- MMPBSA_accuracy.png : a plot of predicted vs. experimental binding affinity results
- MMPBSA_statistics.txt : mean absolute error, root mean square error, and kendall tau of predicted binding affinity results compared to experimentally obtained binding affinities
- rmsd_results.xlsx : ligand RMSD in every AMBER protein-ligand simulation (throughout the shortened trajectory)
- rmsd_box_whisker.png : box-and-whisker plot of ligand RMSDs in every AMBER protein-ligand simulation
