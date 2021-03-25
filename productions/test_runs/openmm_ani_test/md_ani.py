from __future__ import print_function
from simtk.openmm import *
from simtk.openmm.app import *
from parmed.openmm import *
from simtk import unit
from sys import stdout
import sys
from openmmml import MLPotential

prmtop_file=sys.argv[1]
crd_file=sys.argv[2]
prmtop = AmberPrmtopFile(prmtop_file)
inpcrd = AmberInpcrdFile(crd_file)

#potential = MLPotential('ani2x')
#ml_atoms  = [atom.index for atom in prmtop.topology.atoms() if atom.residue.name == "MOL"]
#mm_system = prmtop.createSystem(nonbondedMethod=NoCutoff)
#ml_system = potential.createMixedSystem(prmtop.topology, mm_system, ml_atoms)


potential = MLPotential('ani2x')
ml_atoms  = [atom.index for atom in prmtop.topology.atoms() if atom.residue.name == "MOL289"]
mm_system = prmtop.createSystem(nonbondedMethod=PME,
nonbondedCutoff=1.0*unit.nanometers, constraints=HBonds, rigidWater=True,
ewaldErrorTolerance=0.0005)
ml_system = potential.createMixedSystem(prmtop.topology, mm_system, ml_atoms)

integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)
platform = Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, ml_system, integrator, platform, properties)
# Set the current positions
simulation.context.setPositions(inpcrd.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

simulation.reporters.append(app.DCDReporter('equilibration.dcd', 100))
simulation.reporters.append(app.StateDataReporter(stdout, 100, step = True, potentialEnergy = True, kineticEnergy=True, temperature = True, density = True, volume = True , totalEnergy= True, separator='\t'))

print('Equilibrating...')
simulation.step(100)
