from pyfodmc import pyfodmc
import fodmc 

# Simple test for molecule
sys = 'SO2.xyz'
con_mat = ['(1-2)-(2-2)','(1-3)-(2-2)\n']
pyfodmc.write_pyfodmc_molecules(sys=sys,con_mat=con_mat)
# Fortran call
fodmc.fodmc.get_guess()

