from fodMC.pyfodmc import pyfodmc

# Simple test for molecule
sys = 'SO2.xyz'
con_mat = ['(1-2)-(2-2)','(1-3)-(2-2)\n']
pyfodmc.write_pyfodmc_molecules(sys=sys,con_mat=con_mat)
# Fortran call
output_mode = ['NRLMOL','PyFLOSIC'][1]
output_name = ['',      'SO2_FODs.xyz'][1]
pyfodmc.get_guess(output_mode,output_name)
