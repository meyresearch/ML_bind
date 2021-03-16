from __future__ import print_function
import mdtraj as md
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

'''
Calculate the RMSD for the ligand using a shortened production as input (all 300 frames). RMSD (nm) is plotted vs. Frame (10 ps/frame). Also saves RMSD to a csv file.
This uses "production_short.mdcrd" and "solvated_complex.prmtop" as input, so make sure it is used within the directory of the protein-ligand trajectory that you want to carry out RMSD analysis for.
Usage: python rmsd_analysis.py
'''
#load trajectory
print("Loading trajectory...")
traj_raw = md.load_mdcrd('production_short.mdcrd', top='solvated_complex.prmtop', stride=1)
#remove solvent
traj = traj_raw.atom_slice(traj_raw.top.select('protein or resname MOL'))
#select ligand atoms
lig= traj.top.select('resname MOL')
#calculate rmsd
rmsd = md.rmsd(traj, traj,atom_indices=lig)

#plot results
plt.figure()
plt.plot(rmsd)
plt.title('RMSD to first frame')
plt.xlabel('Frame (10 ps/frame)')
plt.ylabel('RMSD (nm)')
plt.savefig('rmsd.png', dpi=300)
plt.show()
plt.close()

Print("RMSD figure saved as rmsd.png")

rmsd_results = pd.DataFrame(columns=["Frame (10 ps/frame)"])
for i in range(600):
    rmsd_results=rmsd_results.append({"Frame (10 ps/frame)":i}, ignore_index=True)
rmsd_results["RMSD (nm)"]=pd.Series(rmsd)
rmsd_results.to_csv("rmsd_results.csv")

Print("RMSD results saved as rmsd_results.csv")
