import sys
import os
import math
import numpy


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
old = open(pdb_file, 'r')
new = open(out_pdb, 'w')
line = old.readline()
atom_list = []
old_index = 0
max_x = 0
max_y = 0
max_z = 0
index = 0
while line:
    line_split = line.split()
    if len(line_split) > 3:
        if (line_split[0] == 'ATOM' or line_split[0] == 'HETATM') and (
                    line[12:16] == '  CA' or line[12:16] == 'CA  ' or line[12:16] == ' CA '):
            index = line[22:26]
            old_index = line[6:11]
            line2 = line.replace(old_index, index)
            index_int=int(index)
            # check if indicies are in order
            if len(atom_list)>0:
                old_int=atom_list[len(atom_list)-1][0]
                if old_int==index_int-1:
                    pass
                else:
                    print('Warning!!! Went from index '+str(old_int)+' to index '+str(index_int)+'. This is likely a sign that there are multiple chains in the system or a particular atom is not being recognized. Remove multiple chains or excess atoms.')

            type=line[17:20]
            atom_list.append([index_int,type])
            new.write(line2)
            #print(atom_list)

    line = old.readline()
old.close()
new.write('END\n')
new.close()

N=len(atom_list)
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
            for i in range(0,N):
                index=i+1
                new2.write('\t'+str(index)+'\t'+atom_list[i][1]+'\t'+str(index)+'\t'+atom_list[i][1]+'\t'+'CA'+'\t'+str(index)+'\n')

        elif line_split[1] == 'bonds':
            new2.write(line)
            line = template.readline()
            new2.write(line)
            for i in range(0,N-1):
                index=i+1
                next_i=i+2
                new2.write(str(index)+'\t'+str(next_i)+'\t'+'1'+'\t'+str(bond_l)+'\t\t'+str(bond_k)+'\n')

        else:
            new2.write(line)
    else:
        new2.write(line)
    line=template.readline()

new2.close()
template.close()