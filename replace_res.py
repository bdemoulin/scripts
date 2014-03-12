#### Preparation of cobram.parm file, and modification of pdb

import sys

from os.path import basename
from os.path import splitext
from optparse import OptionParser

#### Command line options

parser = OptionParser()

parser.add_option("-s","--solvent",help="Choose the solvent molecule",default="WAT",
                  action="store", dest="solvent")

parser.add_option("-d","--distance",help="Distance from QM part for solvent Medium Layer",default="0.0",
                  action="store", dest="distance")

parser.add_option("-m","--method",help="center of mass of high layer (CoMH), center of mass of high and low layers (CoMHM), within distance to HL atoms (DtH) or within distance to HL/ML atoms (DtHM)",default=" ",
                  action="store", dest="method",choices=("CoMH","CoMHM","DtH","DtHM"," "))

(options, args) = parser.parse_args()

solvent = options.solvent
distance = options.distance
method = options.method

#### Read the initial pdb file and create a modified one
pdb = open(sys.argv[1],"r").readlines()

base = basename(sys.argv[1])
nameMod = str(splitext(base)[0])+"_mod.pdb"

pdb_mod = open(nameMod,"w")

Link1 = 0
Link2 = 0
EndRes = 0

for line in pdb:
    # Modify the first and the last residue 
    if line[17:20] == "THR" and line[22:26] == "   1":
        pdb_mod.write("%sTH1%s" % (line[0:17],line[20:]))
    elif line[17:20] == "LYS" and line[22:26] == " 133":
        pdb_mod.write("%sLY1%s" % (line[0:17],line[20:]))
    # Find the first (Link1) and last (EndRes) atoms of QM part, and the break (Link2)
    elif line[17:20] == "RET" and line[12:16] == " CE ":
        pdb_mod.write("%s" % (line[0:]))
        Link1 = str(line[7:11])
    elif line[17:20] == "RET" and line[12:16] == " CD ":
        pdb_mod.write("%s" % (line[0:]))
        Link2 = str(line[7:11])
    elif line[17:20] == "RET" and line[12:16] == " H39":
        pdb_mod.write("%s" % (line[0:]))
        EndRes = str(line[7:11])
    else:
        pdb_mod.write("%s" % (line[0:]))


#### Create cobram.parm
cobramFile = open("cobram.parm","w")

cobramFile.write("""!!!!!!!!!!!!pdb file 
coordinate file = %s 
 
!!!!!!!!!!!!must be in $AMBERHOME/dat/leap/cmd/ or in a sub-directory therein 
force fields = leaprc.ff99SB_mod  
 
!!!!!!!!!!!!must be in $AMBERHOME/dat/leap/parm/ (check with tleap -f 'force field' which parameter file is used) 
force field parameters = parm99.dat  
 
!!!!!!!!!!!!ATTENTION: 3 inputs needed for each non-standard residue: 1. a 3-letter code 2. structure file (.in, .mol2, .pdb AMBERHOME/dat/leap/prep/); 3. parm file ($AMBERHOME/dat/leap/parm/) 
pre-defined non-standard residues =  
 
!!!!!!!!!!!!3-letter code for non-standard residue to be created (leave empty when standard residues are used) 
non-standard residues =   
        
!!!!!!!!!!!!for the automatic creation of non-standard residues tell antechamber to use either amber or gaff for each residue (leave empty when standard residues are used) 
non-standard residue parameters =  
        
!!!!!!!!!!!!non-stdandard residue charge (leave empty when standard residues are used) 
non-standard residue charges =  
        
!!!!!!!!!!!!3-letter code for solvent (selection according to H-bonding); '0.0' for excluding solvent molecules from medium layer; any number (2.0,3.5, etc.) for all solvent molecules within a radius; center of mass of high layer (CoMH), center of mass of high and low layers (CoMHM), within distance to HL atoms (DtH) or within distance to HL/ML atoms (DtHM) 
solvent = %s %s %s
 
!!!!!!!!!!!!numbering of the high layer atoms as in the pdb file 
high layer atoms = %s-%s  
 
!!!!!!!!!!!!QM charge 
high layer charge = +1 
 
!!!!!!!!!!!!non-solvent medium layer atoms (leave empty in case of a H-L calculation) 
medium layer atoms =  
 
!!!!!!!!!!!!numbering of atom pairs (leave empty if there aren't any) 
links = %s-%s 
 
!!!!!!!!!!!!force additional bonds (exmpl. S-S bonds); naming should be as in tleap (RESID_NUMBER.ATOM_TYPE) (leave empty if there aren't any) 
additional bonds =   
""" % (nameMod, solvent, distance, method, Link1, EndRes, Link1, Link2))
         



   
