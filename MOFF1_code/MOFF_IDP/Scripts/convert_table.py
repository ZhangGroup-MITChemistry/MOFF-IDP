import sys
import os
import math
import numpy

# Potential of the form U=q_i*q_j*exp(-kappa*r)/r+A / r^12-B 0.5*(1+tanh(eta * (rmax-x) ) )
# Expects f(r), -f'(r), g(r), -g'(r), h(r), -h'(r).
# f(r) is charge, g(r) is attractive term, h(r) is repulsive term
# dr=0.002 nm in most tables

# files to export. This script DOES NOT calculate electrostatics, look to add in the future
prev_mat_file=sys.argv[1]
prev_mat=numpy.loadtxt(prev_mat_file)

new_mat_file=sys.argv[2]


# Import R and electostatics to new script
N=numpy.size(prev_mat,0)

new_mat=numpy.zeros((N,7))

new_mat[:,0]=prev_mat[:,0]
new_mat[:,1]=prev_mat[:,1]
new_mat[:,2]=prev_mat[:,2]

skip=8
# parameters for attractive well in nm
eta=7
r_cut=0.8

for i in range(skip,N):
    r=new_mat[i,0]

    g= -0.5 * (1.0 + numpy.tanh( eta*(r_cut-r) ) )
    g2= -0.5 * eta / numpy.power((numpy.cosh( eta*(r_cut-r)  )),2)
    h= numpy.power((1.0/r),12)
    h2=12*numpy.power((1.0/r),13)

    new_mat[i, 3]=g
    new_mat[i,4]=g2
    new_mat[i,5]=h
    new_mat[i,6]=h2



# save new matrix to file
numpy.savetxt(new_mat_file,new_mat)
