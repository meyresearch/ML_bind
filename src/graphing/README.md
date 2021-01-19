# Using 'graphing.py' script:

python graphing.py (path to csv file) (variable you want plotted)

### This script has been produced to plot results in 'equilibration.csv' and 'production.csv' output files from the md.py OpenMM script

### The following strings can be passed as the second argument:
    'totalenergy': produces a graph of total energy (x10^3 kcal/mol) vs. time (ps)
    'kineticenergy': produces a graph of kinetic energy (x10^3 kcal/mol) vs. time (ps)
    'potentialenergy': produces a graph of potential energy (x10^3 kcal/mol) vs. time (ps)
    'temperature': produces a graph of temperature (K) vs. time (ps)
    'density': produces a graph of density (g/(item*mL)) vs. time (ps)
    'volume': produces a graph of volume (nm^3) vs. time (ps)
    'all': produces all the aforementioned graphs

