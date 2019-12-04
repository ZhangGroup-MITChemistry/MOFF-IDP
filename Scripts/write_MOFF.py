import sys
import os
import math
import numpy
import MDAnalysis as mda


if len(sys.argv)>3 or len(sys.argv)==1:
    print("Incorrect usage. Please input proper files:")
    print("python write_ff4.py jobid template_file")
    print("note that the pdb should be the same as the jobid")
    print("Quiting")
    exit()


struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
        pdb = struct_id[-4:]
else:
        pdb = struct_id

# new file types
pdb_file=pdb+'.pdb'
out_pdb=pdb+'_CA.pdb'
top_file=pdb+'.top'

if len(sys.argv)==3:
    template_file=sys.argv[2]
    if template_file[-4:].lower()==".top":
        pass
    else:
        template_file = template_file+".top"
else:
    template_file = 'template_MOFF.top'

# write pdb with only CA
u = mda.Universe(pdb_file,multiframe='no')
calphas = u.select_atoms("name CA")
calphas.write(out_pdb)
# Load new pdb to fix atom indexing
new_universe = mda.Universe(out_pdb,multiframe='no')
calphas= new_universe.select_atoms("name CA")

# debugging
#print(calphas.atoms.n_segments)
#print(calphas.atoms[1].resname)

# Constants
N=len(calphas)
bond_l_i=3.8 # A
bond_k_i=10 # kJ/A^2

# convert units for gromacs
bond_l=0.1*bond_l_i # in nm
bond_k=100*bond_k_i # in kJ/nm^2

# write topology file
new2=open(top_file,'w')
template=open(template_file,'r')
line=template.readline()
while line:
    line_split=line.split()
    if len(line_split)>2:
        if line_split[1] == 'atoms':
            new2.write(line)
            line = template.readline()
            new2.write(line)
            # write all CA atoms
            for i in range(0,N):
                index=calphas.atoms[i].index+1
                res=calphas.atoms[i].resname
                new2.write('\t'+str(index)+'\t'+res+'\t'+str(index)+'\t'+res+'\t'+'CA'+'\t'+str(index)+'\n')

        elif line_split[1] == 'bonds':
            new2.write(line)
            line = template.readline()
            new2.write(line)
            # Write all bonds
            for i in range(0,N-1):
                if calphas.atoms[i].segid==calphas.atoms[i+1].segid:
                    index=calphas.atoms[i].index+1
                    next_i=calphas.atoms[i+1].index+1
                    new2.write(str(index)+'\t'+str(next_i)+'\t'+'1'+'\t'+str(bond_l)+'\t\t'+str(bond_k)+'\n')

        else:
            new2.write(line)
    else:
        new2.write(line)
    line=template.readline()

new2.close()
template.close()