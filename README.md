Scripts
======

Some useful scripts for quantum chemistry or molecular dynamics analysis.

## Description of the scripts

**find_residues**: from a netCDF trajectory and an AMBER topology, finds the number of
residues which are at a given distance of another residue. For instance, it gives for 
each frame the number of water molecules around an organic chromophore. The result is a
normalized histogram. Uses [MDAnalysis toolkit](http://code.google.com/p/mdanalysis/) for 
the analysis and [Prettyplotlib](http://olgabot.github.io/prettyplotlib/) for the plots.

**find_energies.py**: from a MOLCAS *caspt2.output* file, extracts the excitation energies
at RASSCF, CASPT2, and MS-CASPT2 levels. It returns a *csv* table, with the energies in eV 
and in nm.

**replace_res.py**: from a pdb, produces a *cobram.parm* file for COBRAMM computations.

