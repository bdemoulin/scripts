import os
import sys
import csv
import numpy as np

file_to_open = "caspt2.output"

####################################################

def excitation_ev(excited,ground):
    excited = float(excited)
    ground = float(ground)
    excitation_h = excited - ground 
    excitation_ev = 27.2105*excitation_h
    return excitation_ev

####################################################

caspt2_output = open(file_to_open,"r").readlines()

i_root = 0

# Initialization of the tables for the raw energy
RASSCF_energy = []
CASPT2_energy = []
MSCASPT2_energy = []

# Initialization of the tables for the computed excitation energies
RASSCF_excitation = []
CASPT2_excitation = []
MSCASPT2_excitation = []

# Read the file to find the corresponding energies
for line in caspt2_output:
    if str(line).startswith("::"): 
         line = line.split()
         # Save RASSCF Data
         if line[1] == "RASSCF":
             RASSCF_energy += [line[8]]
             i_root += 1
         # Save CASPT2 Data
         if line[1] == "CASPT2":
             CASPT2_energy += [line[6]]
         # Save MS-CASPT2 Data
         if line[1] == "MS-CASPT2":
             MSCASPT2_energy += [line[6]]

# Calculation of excitation energies for RASSCF
for j in range(1, i_root):
    excite_RASSCF = "%.5f" % excitation_ev(RASSCF_energy[j],
                                           RASSCF_energy[0])
    RASSCF_excitation += [ excite_RASSCF ]
    excite_CASPT2 = "%.5f" % excitation_ev(CASPT2_energy[j],
                                           CASPT2_energy[0])
    CASPT2_excitation += [ excite_CASPT2 ]
    excite_MSCASPT2 = "%.5f" % excitation_ev(MSCASPT2_energy[j],
                                             MSCASPT2_energy[0])
    MSCASPT2_excitation += [ excite_MSCASPT2 ]

excitation = np.array( [ (RASSCF_excitation), 
                         (CASPT2_excitation), 
                         (MSCASPT2_excitation) ], 
                         dtype=float )

excitation_nm = 1240 / excitation

# Write excitations in output file
np.savetxt("excitations_nm.csv", excitation_nm, fmt='%.2f')
np.savetxt("excitations.csv", excitation, fmt='%.5f')

# Print the result 
result = os.system("cat excitations.csv")
result_nm = os.system("cat excitations_nm.csv")
