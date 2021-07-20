from fodMC.pyfodmc import pyfodmc

pyfodmc.write_pyfodmc_atoms(sys='C')

output_mode = ['NRLMOL','PyFLOSIC'][0]
output_name = ['',      'C.xyz'][0]
pyfodmc.get_guess(output_mode,output_name)

