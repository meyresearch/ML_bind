from __future__ import print_function
from simtk.openmm import * 
from simtk.openmm.app import *
from parmed.openmm import *
from simtk import unit
from sys import stdout
import sys

def equilibration(prmtop_file, crd_file):
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(crd_file)
    system = prmtop.createSystem(nonbondedMethod=PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=HBonds, rigidWater=True, 
    ewaldErrorTolerance=0.0005)
    integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)
    platform = Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'mixed', 'OpenCLDeviceIndex': str(deviceindex)}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    # Set the current positions
    simulation.context.setPositions(inpcrd.positions)

    # Minimization
    print('Minimizing...')
    simulation.minimizeEnergy()

    # Generate some velocities
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    print('Equilibrating...')

    # Simulation reporters
    simulation.reporters.append(DCDReporter('equilibration.dcd', 100))
    simulation.reporters.append(StateDataReporter('equilibration.csv', 100, step = True, potentialEnergy = True, kineticEnergy=True, temperature = True, density = True, volume = True , totalEnergy= True, separator='\t'))
    simulation.step(5000)

    # Saving data
    positions = simulation.context.getState(getPositions=True).getPositions()
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    PDBFile.writeFile(simulation.topology, positions, open('equilibration.pdb', 'w'))
    checkpoint = 'equilibration.chk'
    simulation.saveCheckpoint(checkpoint)
    return simulation

def production(simulation):

    system = simulation.system
    # NPT for production:
    system.addForce(MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))
    
    #Removing equilibration reporters
    simulation.reporters.pop(-1)
    simulation.reporters.pop(-1)
    
    print("Running production")
    # Simulation reporters
    simulation.reporters.append(DCDReporter('production.dcd', 10000))
    simulation.reporters.append(StateDataReporter('production.csv', 10000, step = True, potentialEnergy = True, kineticEnergy=True, temperature = True, density = True,volume=True, totalEnergy= True, separator='\t'))
    simulation.reporters.append(MdcrdReporter('production.mdcrd', 10000))
    simulation.step(2500000)
    
    # Save final frame to PDB file
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open('production.pdb', 'w'))


if __name__ == '__main__':
    deviceindex = sys.argv[3]
    simulation = equilibration(sys.argv[1], sys.argv[2])
    production(simulation)




