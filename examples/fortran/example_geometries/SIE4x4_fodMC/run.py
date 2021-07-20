import fodmc
# output_mode: PyFLOSIC, NRLMOL
# output_name: NameOfMolecule.xyz  (for PyFLOSIC only)
output_mode = ['NRLMOL','PyFLOSIC'][1]
output_name = ['',      'Nuc_FOD.xyz'][1]
fodmc.fodmc_mod.get_guess(output_mode,output_name)
