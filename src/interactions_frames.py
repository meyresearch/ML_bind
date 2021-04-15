from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import seaborn as sns
import argparse
import sys

'''
Description: Calculate the percentage of frames in which certain common protein-ligand interactions occur,
e.g. Van Der Waals' contacts, hydrogen bonding, and a cation-pi interaction for ejm_49.
Results are plotted as a graph (interactions_frames.png) and saved as a .csv file (interactions_frames.csv)
If the script is being run on the complex using ejm_49, please use the --ejm_49 option to calculate
cation-pi interactions. This script must be used within a production simulation folder, as it reads the
Usage: python atomic_distances.py [--ejm_49]
'''

parser = argparse.ArgumentParser()
parser.add_argument("--ejm_49", help="add cation-pi interactions for ligand ejm_49", action="store_true")
args = parser.parse_args()

#Load trajectory
print("Loading trajectory...")
traj = md.load_mdcrd('production_short.mdcrd', top='solvated_complex.prmtop', stride=1)

print("Calculating atomic distances...")
#Note: See project paper for reasoning of cutoff distance used for various interactions

#Common Van Der Waals' Contacts:
#Cutoff of 4.0 angstroms used for VDW contacts

#Calculate distance between Val92 atom CG2 (sidechain C) and ligand atom C10 (pyridine) in each frame
vdw_pyridine1_val92=traj.top.select("resid 91 and name CG2")
vdw_pyridine1_mol289=traj.top.select("resid 288 and name C10")
vdw_pyridine1_indices=np.array([vdw_pyridine1_val92[0], vdw_pyridine1_mol289[0]], ndmin=2)
vdw_pyridine1=md.compute_distances(traj, vdw_pyridine1_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_pyridine1=vdw_pyridine1[vdw_pyridine1<=0.4]

#Calculate distance between Val22 atom CG1 (sidechain C) and ligand atom C5 (dichlorobenzene) in each frame
vdw_dichlorobenzene1_val22=traj.top.select("resid 21 and name CG1")
vdw_dichlorobenzene1_mol289=traj.top.select("resid 288 and name C5")
vdw_dichlorobenzene1_indices=np.array([vdw_dichlorobenzene1_val22[0], vdw_dichlorobenzene1_mol289[0]], ndmin=2)
vdw_dichlorobenzene1=md.compute_distances(traj, vdw_dichlorobenzene1_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_dichlorobenzene1=vdw_dichlorobenzene1[vdw_dichlorobenzene1<=0.4]

#Calculate distance between Ala39 atom CB (sidechain C) and ligand atom C10 (pyridine) in each frame
vdw_pyridine2_ala39=traj.top.select("resid 38 and name CB")
vdw_pyridine2_mol289=traj.top.select("resid 288 and name C10")
vdw_pyridine2_indices=np.array([vdw_pyridine2_ala39[0], vdw_pyridine2_mol289[0]], ndmin=2)
vdw_pyridine2=md.compute_distances(traj, vdw_pyridine2_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_pyridine2=vdw_pyridine2[vdw_pyridine2<=0.4]

#Calculate distance between Ile71 atom CD1 (sidechain C) and ligand atom C9 (pyridine) in each frame
vdw_pyridine3_ile71=traj.top.select("resid 70 and name CD1")
vdw_pyridine3_mol289=traj.top.select("resid 288 and name C9")
vdw_pyridine3_indices=np.array([vdw_pyridine3_ile71[0], vdw_pyridine3_mol289[0]], ndmin=2)
vdw_pyridine3=md.compute_distances(traj, vdw_pyridine3_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_pyridine3=vdw_pyridine3[vdw_pyridine3<=0.4]

#Calculate distance between Leu141 atom CD2 (sidechain C) and ligand atom C3 (dichlorobenzene) in each frame
vdw_dichlorobenzene2_leu141=traj.top.select("resid 140 and name CD2")
vdw_dichlorobenzene2_mol289=traj.top.select("resid 288 and name C3")
vdw_dichlorobenzene2_indices=np.array([vdw_dichlorobenzene2_leu141[0], vdw_dichlorobenzene2_mol289[0]], ndmin=2)
vdw_dichlorobenzene2=md.compute_distances(traj, vdw_dichlorobenzene2_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_dichlorobenzene2=vdw_dichlorobenzene2[vdw_dichlorobenzene2<=0.4]

#Calculate distance between Leu141 atom CD2 (sidechain C) and ligand atom C4 (dichlorobenzene) in each frame
vdw_dichlorobenzene3_leu141=traj.top.select("resid 140 and name CD2")
vdw_dichlorobenzene3_mol289=traj.top.select("resid 288 and name C4")
vdw_dichlorobenzene3_indices=np.array([vdw_dichlorobenzene3_leu141[0], vdw_dichlorobenzene3_mol289[0]], ndmin=2)
vdw_dichlorobenzene3=md.compute_distances(traj, vdw_dichlorobenzene3_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_dichlorobenzene3=vdw_dichlorobenzene3[vdw_dichlorobenzene3<=0.4]

#Calculate distance between Leu14 atom CD2 (sidechain C) and ligand atom C13 (pyridine) in each frame
vdw_pyridine4_leu14=traj.top.select("resid 13 and name CD1")
vdw_pyridine4_mol289=traj.top.select("resid 288 and name C11")
vdw_pyridine4_indices=np.array([vdw_pyridine4_leu14[0], vdw_pyridine4_mol289[0]], ndmin=2)
vdw_pyridine4=md.compute_distances(traj, vdw_pyridine4_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_pyridine4=vdw_pyridine4[vdw_pyridine4<=0.4]

#Calculate distance between Asp152 atom CD2 (sidechain C) and ligand atom C6 (dichlorobenzene) in each frame
vdw_dichlorobenzene4_asp152=traj.top.select("resid 151 and name CB")
vdw_dichlorobenzene4_mol289=traj.top.select("resid 288 and name C6")
vdw_dichlorobenzene4_indices=np.array([vdw_dichlorobenzene4_asp152[0], vdw_dichlorobenzene4_mol289[0]], ndmin=2)
vdw_dichlorobenzene4=md.compute_distances(traj, vdw_dichlorobenzene4_indices)
#Apply cutoff and only save distances within cutoff range
cutoff_vdw_dichlorobenzene4=vdw_dichlorobenzene4[vdw_dichlorobenzene4<=0.4]


#Common Hydrogen bonds:

#Calculate distance between Val92 atom N (backbone N) and ligand atom N2 (amide N) in each frame
hbond_pyridineN_val92=traj.top.select("resid 91 and name N")
hbond_pyridineN_mol289=traj.top.select("resid 288 and name N2")
hbond_pyridineN_indices=np.array([hbond_pyridineN_val92[0], hbond_pyridineN_mol289[0]], ndmin=2)
hbond_pyridineN=md.compute_distances(traj, hbond_pyridineN_indices)
#Apply cutoff of 3.23 angstroms for heavy atoms in N-HN hydrogen bonds and only save distances within cutoff range
cutoff_hbond_pyridineN=hbond_pyridineN[hbond_pyridineN<=0.323]

#Calculate distance between Val92 atom O (backbone O) and ligand atom N3 (pyridine N) in each frame
hbond_amideH_val92=traj.top.select("resid 91 and name O")
hbond_amideH_mol289=traj.top.select("resid 288 and name N3")
hbond_amideH_indices=np.array([hbond_amideH_val92[0], hbond_amideH_mol289[0]], ndmin=2)
hbond_amideH=md.compute_distances(traj, hbond_amideH_indices)
#Apply cutoff of 3.17 angstroms for heavy atoms in O-HN hydrogen bonds and only save distances within cutoff range
cutoff_hbond_amideH=hbond_amideH[hbond_amideH<=0.317]


#Cation-pi interaction for Arg12 with ligand ejm_49:


if args.ejm_49:
    #Calculate distance between Arg12 atom CZ (carbamimidamido carbon) and ligand atom C14 (benzyl carbon)
    cation_pi_arg12= traj.top.select("resid 11 and name CZ")
    cation_pi_mol289= traj.top.select("resid 288 and name C16")
    cation_pi_indices= np.array([cation_pi_arg12[0], cation_pi_mol289[0]], ndmin=2)
    cation_pi = md.compute_distances(traj, cation_pi_indices)
#Apply cutoff distance of 5 angstroms for cation-pi interactions
    cutoff_cation_pi=cation_pi[cation_pi<=0.5]


#Create dictionary of results:
    results={'Interaction':["Val92(CG2)-Lig(C10)", "Val22(CG1)-Lig(C5)",
                        "Ala39(CB)-Lig(C10)", "Ile71(CD1)-Lig(C9)",
                        "Leu141(CD2)-Lig(C3)","Leu141(CD2)-Lig(C4)",
                        "Leu14(CD1)-Lig(C11)", "Asp152(CB)-Lig(C6)",
                        "Val92(NH)-Lig(Pyridine N)",
                        "Val92(O)-Lig(Amide NH)", "Arg12(Cz)-Lig(Phenyl C16)"],
         'Frames Present (%)':[np.shape(cutoff_vdw_pyridine1)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene1)[0]/6,
                           np.shape(cutoff_vdw_pyridine2)[0]/6,
                           np.shape(cutoff_vdw_pyridine3)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene2)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene3)[0]/6,
                           np.shape(cutoff_vdw_pyridine4)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene4)[0]/6,
                           np.shape(cutoff_hbond_pyridineN)[0]/6,
                           np.shape(cutoff_hbond_amideH)[0]/6,
                           np.shape(cutoff_cation_pi)[0]/6]}

else:
    results={'Interaction':["Val92(CG2)-Lig(C10)", "Val22(CG1)-Lig(C5)",
                        "Ala39(CB)-Lig(C10)", "Ile71(CD1)-Lig(C9)",
                        "Leu141(CD2)-Lig(C3)","Leu141(CD2)-Lig(C4)",
                        "Leu14(CD1)-Lig(C11)", "Asp152(CB)-Lig(C6)",
                        "Val92(NH)-Lig(Pyridine N)",
                        "Val92(O)-Lig(Amide NH)"],
         'Frames Present (%)':[np.shape(cutoff_vdw_pyridine1)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene1)[0]/6,
                           np.shape(cutoff_vdw_pyridine2)[0]/6,
                           np.shape(cutoff_vdw_pyridine3)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene2)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene3)[0]/6,
                           np.shape(cutoff_vdw_pyridine4)[0]/6,
                           np.shape(cutoff_vdw_dichlorobenzene4)[0]/6,
                           np.shape(cutoff_hbond_pyridineN)[0]/6,
                           np.shape(cutoff_hbond_amideH)[0]/6]}

#Create dataframe of results:
df_results=pd.DataFrame(results)
df_results.to_csv("interactions_frames.csv", index=False)
print("Results saved as 'interactions_frames.csv'")

#Plot results to bar chart
sns.set_theme(context='paper', style='white', font_scale=2)
sns.set_style("ticks")
colours=["green", "green", "green", "green", "green", "green", "green", "green", "blue", "blue", "red"]
df_results.plot.bar(x="Interaction",y="Frames Present (%)", color=colours, figsize=(20,10), fontsize=25, legend=None)
plt.xlabel("Interaction", fontsize=30)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Frames Present (%)", fontsize=30)
sns.despine()
plt.tight_layout()
#Save plot
plt.savefig("interactions_frames.png", dpi=300)

print("Graph saved as 'interactions_frames.png'")
