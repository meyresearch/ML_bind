import BioSimSpace as BSS

from Sire.Mol import AtomIdx

import os
import sys

# Extract the directory for the job.
try:
    job_dir = os.getenv("JOB_DIR")
except:
    job_dir = None

# No job directory set, use the current directory.
if job_dir is None:
    job_dir = "."

# Extract the ligand numbers.
num0 = sys.argv[1]
num1 = sys.argv[2]

# Load the protein and crystal waters.
protein = BSS.IO.readMolecules(BSS.IO.glob("%s/protein/TYK2/parameters/vacuum*" % job_dir)).getMolecules()[0]

# Extract the waters.
# waters = protein_water.getWaterMolecules()

# Parameterise the protein.
# protein = BSS.Parameters.ff14SB(protein_water.getMolecules()[0]).getMolecule()

# Load the parameterised ligands.
lig0 = BSS.IO.readMolecules(BSS.IO.glob("%s/poses/parametrised/tyk_lig%s/vacuum.*" % (job_dir, num0))).getMolecules()[0]
lig1 = BSS.IO.readMolecules(BSS.IO.glob("%s/poses/parametrised/tyk_lig%s/vacuum.*" % (job_dir, num1))).getMolecules()[0]

# If a mapping file exists, then load the mapping. Otherwise, use BioSimSpace
# to create the mapping.
mapping = {}

# Forward mapping.
if os.path.isfile("%s/FESetup_mappings/merge_errors/%s_%s.txt" % (job_dir, num0, num1)):
    with open("%s/FESetup_mappings/merge_errors/%s_%s.txt" % (job_dir, num0, num1), "r") as file:
        for line in file:
            pair = line.strip().split()
            mapping[AtomIdx(int(pair[0]))] = AtomIdx(int(pair[1]))

# Reverse mapping.
elif os.path.isfile("%s/FESetup_mappings/merge_errors/%s_%s.txt" % (job_dir, num1, num0)):
    with open("%s/FESetup_mappings/merge_errors/%s_%s.txt" % (job_dir, num1, num0), "r") as file:
        for line in file:
            pair = line.strip().split()
            # Invert the indices.
            mapping[AtomIdx(int(pair[1]))] = AtomIdx(int(pair[0]))

# No mapping, generate it ourselves.
else:
    # Find the best mapping of atoms between the ligands.
    mapping = BSS.Align.matchAtoms(lig0, lig1)

# Align lig0 to lig1 based on the mapping.
lig0 = BSS.Align.rmsdAlign(lig0, lig1, mapping)

# Merge the two ligands based on the mapping.
merged = BSS.Align.merge(lig0, lig1, mapping, allow_ring_breaking=True)

# Create the composite system.
system = merged + protein

# Solvate in a 60 angstrom box of TIP3P water.
solvated = BSS.Solvent.tip3p(molecule=system, box=3*[80*BSS.Units.Length.angstrom])

# Create the free energy protocol.
protocol = BSS.Protocol.FreeEnergy(timestep=2*BSS.Units.Time.femtosecond, runtime=4*BSS.Units.Time.nanosecond, num_lam=12)

# Initialise the binding free energy object.
freenrg = BSS.FreeEnergy.Binding(solvated, protocol, work_dir="TYK2_%s_%s" % (num0, num1), engine="SOMD", box=3*[30*BSS.Units.Length.angstrom])

# Run the simulation.
freenrg.run()
