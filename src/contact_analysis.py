from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import mdtraj.testing
import numpy as np
from contact_map import ContactTrajectory, RollingContactFrequency, ContactFrequency, ContactDifference

'''
Find hydrogen bonds and contacts between ligand and protein. Inputs are trajectory (mdcrd) file and topology (top, prmtop or prm7) file.

Usage: python contact_analysis.py trajectory_short.mdcrd solvated_complex.top
'''

def contact_analysis(mdcrd_file, top_file):

    #Load topology and trajectory
    topology = md.load_prmtop(top_file)
    print("Loading trajectory...")
    traj = md.load_mdcrd(mdcrd_file, top= topology)

    #Analyse hydrogen bonds between ligand (residue MOL289) and protein
    print("Analysing Hbonds...")
    hbonds = md.baker_hubbard(traj, periodic=False)
    label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
    substring="MOL289"

    print(" ") #Add empty line for clarity
    print("Hydrogen bonds between ligand and protein:")
    print( "   (If empty there were no Hbonds to report)")
    for hbond in hbonds:
        if substring in label(hbond):
            print(label(hbond))
    print(" ")

    #Analyse atom contacts between ligand (residue 288) and protein
    print("Analysing contacts...")
    lig_top = topology.residue(288)
    trajectory_contacts = ContactFrequency(traj)

    print(" ")
    print("Most common residue contacts between residue and protein:")
    print("   (If empty there were no contacts to report)")
    for residue_contact in trajectory_contacts.residue_contacts.most_common(lig_top):
        if residue_contact[1] > 0.8:
            print(residue_contact)
    print(" ")

    print("Most common atom contacts between residue and protein:")
    print("   (If empty there were no contacts to report)")
    for atom_contact in trajectory_contacts.most_common_atoms_for_residue(lig_top):
        if atom_contact[1] > 0.9:
            print(atom_contact)
    print(" ")

if __name__ == '__main__':
    contact_analysis(sys.argv[1], sys.argv[2])
