from fodMC.pyfodmc import pyfodmc

# Simple test for molecule
sys = 'SO2.xyz'
con_mat = ['(1-2)-(2-2)','(1-3)-(2-2)\n']
pyfodmc.write_pyfodmc_molecules(sys=sys,con_mat=con_mat)
# Fortran call
pyfodmc.get_guess()

