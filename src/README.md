# src

This directory contains scripts, jupyter notebooks, or example command line inputs for setting up systems, running simulations, and analysing results. Please see each individual script for comments describing how to use it. 

The following files can be found:
- md.py : OpenMM script used for running protein-ligand simulations
- graphing.py : script used for plotting system conditions (e.g. temperature, energies, etc.) over the course of the equilibration period and saving the results as images
- mmpbsa.in : MM/PBSA input file that was passed to the MMPBSA.py script during binding affinity analysis
- complex2img.py : OpenEye Toolkits script for analysing protein-ligand interactions 
- rmsd_analysis.py : script for calculating the ligand RMSD. Plots and saves results as .png and .csv files
- interactions_frames.py : script for determining the percentage of frames in which a specific protein-ligand contact is present. Plots and saves results as .png and .csv files

The 'examples' folder contains example command line inputs (as .txt files) for preparing structure files for AMBER, creating simulation systems, processing trajectories, analysing results and more. Please see each individual file for a description. The folder also contains Jupyter notebooks for plotting MM/PBSA results (in the 'MMPBSA_plotting' subfolder) and creating box-and-whisker plots from ligand RMSD results (in the 'box_whisker_plot' subfolder).
