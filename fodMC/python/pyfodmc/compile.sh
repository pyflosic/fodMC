# Author: S. Schwalbe 
# Date: 29.04.2019 (also the date when converting fodmc_alex.f90 to fodmc_mod.f90) 
# Note:
#       fodmc_mod.f90 is copy of fodmc_alex.f90, where the fodmc program 
#       is replaced by module containing all subroutines. This is done 
#       to enable the fodmc to be compiled with f2py and used as libaray. 
#       (1) replace (header) 'program fodMC' with 
#           'module fodmc
#           contains
#           subroutine points_on_sphere_metropolis_spin_centers()' 
#       (2) replace (footer) 'end program fodMC' (before subroutine mc_step) with 
#           'end subroutine points_on_sphere_metropolis_spin_centers'
#       (3) add at the end of the file 'end module fodmc'
#                  

# convert fortran rountine to a .so object which can be used by python 
f2py -c -m fodmc fodMC_motif.f90
mv fodmc.*.so fodmc.so
