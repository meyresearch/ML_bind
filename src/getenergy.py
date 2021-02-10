from __future__ import print_function
from simtk.openmm import *
from simtk.openmm.app import *
from parmed.openmm import *
from simtk import unit
from sys import stdout
import sys
import mdtraj as md

'''
Definition:
Get potential energy (in kJ/mol) of the ligand in the last frame of an MD simulation. This only works with simulations containing the protein Tyk2 (as the ligand index is residue 288)

Usage:
python getenergy.py production.pdb (ligand_topology.top) production_short.mdcrd
'''

def get_energy(pdb_file, prmtop_file, mdcrd_file):
    
    print('Loading trajectory...')
    #Load MD trajectory
    full_trajectory = md.load_mdcrd(mdcrd_file,top = pdb_file)
    
    print('Parsing trajectory...')
    #Take snapshot of the last frame, then separate ligand snapshot
    last_frame = full_trajectory[-1]
    ligand = last_frame.topology.select('resid 288')
    lig_traj = last_frame.atom_slice(ligand, inplace=False)
    
    #Check that parsed trajectory is correct, with correct number of atoms (Note that the ligand 'residue' will be defined as 'Mol')
    print('The selected trajectory is:')
    print(lig_traj.topology)
    print('All atoms: %s' % [atom for atom in lig_traj.topology.atoms])
    
    #Save ligand trajectory
    lig_traj.save_amberrst7('lig_last_frame.rst7')
    print('Ligand trajectory saved as lig_last_frame.rst7')
    
    #Create system
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile('lig_last_frame.rst7')
    system = prmtop.createSystem(nonbondedMethod=NoCutoff,
    nonbondedCutoff=1.0*unit.nanometers, constraints=HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)
    integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    simulation = Simulation(prmtop.topology, system, integrator)
    
    # Set the current positions
    simulation.context.setPositions(inpcrd.positions)
    
    #Get potential energy
    potentialenergy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print('The potential energy of the ligand is: ', potentialenergy)
    
    return potentialenergy

if __name__=='__main__':
    get_energy(sys.argv[1], sys.argv[2], sys.argv[3])
