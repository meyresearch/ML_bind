##########################################################################
# this script was generated by openmm-builder. to customize it further,
# you can save the file to disk and edit it with your favorite editor.
##########################################################################

from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import sys

prmtop_file=sys.argv[1]
crd_file=sys.argv[2]
deviceindex = sys.argv[3]
prmtop = app.AmberPrmtopFile(prmtop_file)
inpcrd = app.AmberInpcrdFile(crd_file)
system = prmtop.createSystem(nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(0.00001)
platform = mm.Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'mixed', 'OpenCLDeviceIndex': str(deviceindex)}
simulation = app.Simulation(prmtop.topology, system, integrator, platform,
    properties)
simulation.context.setPositions(inpcrd.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

#simulation.reporters.append(app.DCDReporter('equilibration.dcd', 100))
#simulation.reporters.append(app.StateDataReporter('equilibration.csv', 100, step = True, potentialEnergy = True, kineticEnergy=True, temperature = True, density = True, volume = True , totalEnergy= True, separator='\t'))

print('Equilibrating...')
simulation.step(100)

system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))

simulation.reporters.append(app.DCDReporter('trajectory.dcd', 100))
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True, volume=True,
    speed=True, totalSteps=1000, separator='\t'))

print('Running Production...')
simulation.step(1000)
print('Done!')
