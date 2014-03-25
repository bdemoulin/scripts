# Hydration of the retinal chromophore

# Perform the analysis
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.KDTree.NeighborSearch import AtomNeighborSearch
from MDAnalysis.coordinates.TRJ import NCDFReader

# Required for plots
import matplotlib.pyplot as plt
import matplotlib.pylab as plab
import prettyplotlib as ppl

##########################################

# Parameters: We want to find how many "residues_to_find" are around "radius" 
# of "center"
center = "bynum 1745:1797"
residue_to_find = "resname WAT"
radius = 4.0

# Name of topology (topol.prmtop) and trajectory (mdcrd.ncdf)
# !!!!!! ".prmtop" and ".netcdf" extensions required !!
topology = 'topol.prmtop'
intraj   = 'traj.ncdf'

# Title for the file containing the hydration number
file_title = "hydration_file"

#Parameters for the plot
title = "Hydration for all-free M10"
title_of_file = "hydration_all-free_M10.png"


###########################################

def number_of_neighbouring_residues(residue_to_find, center, radius):
    x = AtomNeighborSearch(residue_to_find)
    neighboring = AtomNeighborSearch.search_list(x, center, radius, level='R')
    number = len(neighboring)
    return number

def plot_figure(title, title_of_file):
    fig, ax = plt.subplots(1)
    ax.set_title(title)
    ppl.bar(bin_edges[:-1], 
            hist, 
            width = 1, 
            xticklabels = histogram_bin,
            annotate = True)
    fig.savefig(title_of_file)

##########################################

# Load a universe, i.e. an object that contains the topology and all the
# available coordinates:
u = Universe(topology,intraj)

# Where are the residues we are interested in in the Universe ?
center_in_file = u.selectAtoms(center)
residue_in_file = u.selectAtoms(residue_to_find)

# Now, we build an histogram
#   - hydration_number_histogram_bin  contains the possible values for the
#     hydration number.
#   - hydration_number stores all the values along the trajectory
histogram_bin = []
hydration_number = np.array([])
num_frame = 0

# File where we write the hydration number
hydration_file = open(file_title,"w")

# On each frame, we find the residues around the center.
# Then, append the array containing the number of residues aroud the center.
# If this number has not appeared before, we store it in histogram_bin
for ts in u.trajectory:
    num_frame += 1
    number_of_residues = number_of_neighbouring_residues(residue_in_file, 
                                                       center_in_file, 
                                                       radius)
    hydration_number = np.append(hydration_number,[number_of_residues])
    hydration_file.write("%d  %d\n" % (num_frame,number_of_residues))
    if number_of_residues not in histogram_bin:
        histogram_bin += [number_of_residues]

print hydration_number

# Construction of histogram in Python requires to add an unreached maximum:
histogram_bin += [max(histogram_bin) + 1]

# Next, we sort the bins:
histogram_bin.sort()

# We then construct the histogram. Here, we have a normalized histogram
hist, bin_edges = np.histogram(hydration_number, 
                               bins = histogram_bin, 
                               normed=True)

# Finally, we plot the histogram
plot_figure(title, title_of_file)


