#cp ../src/fodMC_alex.f90 fodmc.f90
# convert fortran rountine to a .so object which can be used by python 
f2py -c -m fodmc fodmc.f90
