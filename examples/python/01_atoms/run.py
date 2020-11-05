from fodMC.pyfodmc import pyfodmc

pyfodmc.write_pyfodmc_atoms(sys='Ne')

output_mode = ['NRLMOL','PyFLOSIC'][1]
output_name = ['',      'Ne.xyz'][1]
pyfodmc.get_guess(output_mode,output_name)

