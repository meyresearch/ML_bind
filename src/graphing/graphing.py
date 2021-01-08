# Produce a graph of temperature (K) vs. timestep (ps)

#Importing necessary modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sys import stdout
import sys

# Read equilibration csv file
def read_file(csv_file):
    data = pd.read_csv(csv_file, sep= '\t')
    return data

equilibration_data = read_file(sys.argv[1])

# Create timestep and time arrays
timestep = equilibration_data['#"Step"'].to_numpy()
time = timestep/500 #in picoseconds

# Add "Time" column to equilibration dataframe
equilibration_data["Time (ps)"] = time

# Create more arrays of variables in columns
potentialenergy = equilibration_data["Potential Energy (kJ/mole)"].to_numpy()
kineticenergy = equilibration_data["Kinetic Energy (kJ/mole)"].to_numpy()
totalenergy = equilibration_data["Total Energy (kJ/mole)"].to_numpy()
temperature = equilibration_data["Temperature (K)"].to_numpy()
density = equilibration_data["Density (g/mL)"].to_numpy()
volume = equilibration_data["Box Volume (nm^3)"].to_numpy()

# Print dataframe
print(equilibration_data)

# Create plots
def create_plots(y_variable):

    if y_variable == 'totalenergy':
    # Create total energy plot. Saved as "totalenergy.png"
        plt.plot(time, totalenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Total Energy (x$10^3$ kJ/mol)")
        plt.savefig('totalenergy.png', dpi=300)
        plt.clf()

    elif y_variable == 'kineticenergy':
    # Create kinetic energy plot. Saved as "kineticenergy.png"
        plt.plot(time, kineticenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Kinetic Energy (x$10^3$ kJ/mol)")
        plt.savefig('kineticenergy.png', dpi=300)
        plt.clf()

    elif y_variable == 'potentialenergy':
    # Create potential energy plot. Saved as "potentialenergy.png"
        plt.plot(time, potentialenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Potential Energy (x$10^3$ kJ/mol)")
        plt.savefig('potentialenergy.png', dpi=300)
        plt.clf()

    elif y_variable == 'temperature':
    # Create temperature plot. Saved as "temperature.png"
        plt.plot(time, temperature, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")
        plt.savefig('temperature.png', dpi=300)
        plt.clf()

    elif y_variable == 'density':
    # Create density plot. Saved as "density.png"
        plt.plot(time, density, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Density (g/mL)")
        plt.savefig('density.png', dpi=300)
        plt.clf()

    elif y_variable == 'volume':
    # Create volume plot. Saved as "volume.png"
        plt.plot(time, volume, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Box Volume ($nm^3$)")
        plt.savefig('volume.png', dpi=300)
        plt.clf()
    
    elif y_variable == 'all':
        # Create total energy plot. Saved as "totalenergy.png"
        plt.plot(time, totalenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Total Energy (x$10^3$ kJ/mol)")
        plt.savefig('totalenergy.png', dpi=300)
        plt.clf()
        
        # Create kinetic energy plot. Saved as "kineticenergy.png"
        plt.plot(time, kineticenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Kinetic Energy (x$10^3$ kJ/mol)")
        plt.savefig('kineticenergy.png', dpi=300)
        plt.clf()
        
        # Create potential energy plot. Saved as "potentialenergy.png"
        plt.plot(time, potentialenergy/1000, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Potential Energy (x$10^3$ kJ/mol)")
        plt.savefig('potentialenergy.png', dpi=300)
        plt.clf()
        
        # Create temperature plot. Saved as "temperature.png"
        plt.plot(time, temperature, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")
        plt.savefig('temperature.png', dpi=300)
        plt.clf()
        
        # Create density plot. Saved as "density.png"
        plt.plot(time, density, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Density (g/mL)")
        plt.savefig('density.png', dpi=300)
        plt.clf()
        
        # Create volume plot. Saved as "volume.png"
        plt.plot(time, volume, color = 'black', linewidth = '1.5')
        plt.xlabel("Time (ps)")
        plt.ylabel("Box Volume ($nm^3$)")
        plt.savefig('volume.png', dpi=300)
        plt.clf()

create_plots(sys.argv[2])
