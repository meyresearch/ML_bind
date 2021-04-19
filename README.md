# ML_bind

This respository contains all the raw data, simulation files, results, and analysis files for the project. Please follow the directories and guidance below to find the files you are looking for. Please read the README.md file within each directory for further guidance.

#### data:
Contains all the protein and ligand structure files used to set up systems for simulations.

#### docs:
Contains results that are applicable to all protein-ligand systems. If you are looking for results of a specific system, please visit the relevant folder in the 'production' directory.

#### productions:
Contains folders corresponding to each protein-ligand simulation. Each folder contains files including:
- Processed structure files for creating the system
- Topology and input coordinate files for running simulations and conducting MM/PBSA analysis
- Graphs of certain variables over the equilibration period of simulations
- Simulation reporters for equilibration and production periods of simulations
- PDB structure file of the final production simulation frame
- MM/PBSA results
- Ligand RMSD analysis results
- OpenEye Toolkits results for protein-ligand interactions
- MDTraj contact analysis results


(Please note that, due to GitHub file size restrictions, the trajectory files from production simulations have not been added. If you require these trajectory files, please contact the authors.)

#### src:
Contains all the scripts used for running simulations and analysing results. Also contains example command line inputs or jupyter notebooks for preparing systems, conducting MM/PBSA, and analysing results.

