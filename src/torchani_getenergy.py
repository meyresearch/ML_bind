import numpy as np
import pandas as pd
import torch
import torchani

'''
Usage:
python get_ani_energy.py ligand.pdb
'''

def get_ani_energy(pdb_file):
    #Read ligand.pdb file
    ligand = pd.read_csv(pdb_file, delim_whitespace=True, skiprows=1, header =None, usecols=[2,6,7,8,11], names=['atom_id', 'x', 'y', 'z', 'atom_type'])
    ligand_clean = ligand.dropna(thresh=4)
    
    #Get positions and species, format correctly for Torch
    positions=ligand_clean[['x','y','z']].to_numpy()
    atom_types=ligand_clean['atom_type'].to_numpy()
    atom_dictionary = {'H':1, 'C':6, 'N':7, 'O':8, 'F':9, 'P':15, 'S':16, 'Cl':17}
    atoms = np.array([atom_dictionary[letter] for letter in atom_types])
    
    #Set up moedlling system
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = torchani.models.ANI2x(periodic_table_index=True).to(device)
    
    #Create tensors for coodinates and species
    coordinates = torch.tensor([positions], requires_grad=True, device=device).float()
    species = torch.tensor([atoms], device=device)
    
    #Calculate energy, forces
    energy = model((species, coordinates)).energies
    derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
    force = -derivative
    
    #Print energy and force results
    print('Energy:', energy.item(), 'Hartrees')
    print('Force:', force.squeeze())
    
if __name__=='__main__':
    get_ani_energy(sys.argv[1])
