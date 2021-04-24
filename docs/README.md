# docs

This directory contains documents with results of the study that are applicable to all protein-ligand simulations.

If you are looking for results that are specific to a certain ligand (e.g. protein-ligand interactions, specific ligand RMSD plot, etc.), please visit the relevant subfolder for the ligand in the 'productions' directory.

This directory contains the following files:
- MMPBSA_results_amber.csv : contains MM/PBSA binding affinity results for each AMBER protein-ligand simulation, as well as experimental binding affinity results obtained from Wang et al.
- MMPBSA_results_ani.csv : contains MM/PBSA binding affinity results for each ANI-2x/AMBER protein-ligand simulation, as well as experimental binding affinity results obtained from Wang et al.
- MMPBSA_accuracy_amber.png : a plot of predicted vs. experimental binding affinity results. Predicted binding affinities were obtained from AMBER simulations.
- MMPBSA_accuracy_ani.png : a plot of predicted vs. experimental binding affinity results. Predicted binding affinities were obtained from ANI-2x/AMBER simulations.
- MMPBSA_statistics.txt : mean absolute error, root mean square error, and kendall tau of predicted binding affinity results (for both AMBER and ANI-2x/AMBER simulations) compared to experimentally obtained binding affinities. 
- rmsd_results.xlsx : ligand RMSD in every AMBER and ANI-2x/AMBER simulation (throughout the shortened trajectory)
- rmsd_box_whisker_amber.png : box-and-whisker plot of ligand RMSDs in every AMBER simulation
- rmsd_box_whisker_ani.png : box-and-whisker plot of ligand RMSDs in every ANI-2x/AMBER simulation and the subset of AMBER simulations for the same protein-ligand systems
