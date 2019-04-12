# pyfodmc - python interface for fodmc code to pyflosic 
# Author:       Sebastian Schwalbe (SS)  
# Changelog:    11.02.2019 
#               SS update to the new fodmc.f90 
#               SS add new write_pyfodmc_molecules function 
import fodmc as fodmc 
from ase.io import read,write  

def write_pyfodmc_atoms(sys,mc_steps=10000,mc_size=0.001):
    #
    # write fodmc input for atoms 
    #
    # create system file 
    f = open('system','w')
    f.write('1 %s\n' % sys)
    f.write('bohr\n')
    f.write('%s 0.0 0.0 0.0\n' %sys)
    f.write('simulation (MC steps, step size)\n')
    f.write('%1.f\n'  %mc_steps)
    f.write('%0.8f\n' %mc_size)
    f.close()

def write_pyfodmc_molecules(sys,con_mat):
    
    ase_atoms = read(sys)
    sym = ase_atoms.get_chemical_symbols()
    pos = ase_atoms.get_positions() 
    natoms = len(ase_atoms.get_chemical_symbols())
    
    # create system file 
    f = open('system','w')
    f.write('%i %s\n' % (natoms,sys))
    f.write('angstrom\n')
    for p in range(len(pos)):
        f.write('%s %0.5f %0.5f %0.5f\n' %(sym[p],pos[p][0],pos[p][1],pos[p][2]))
    f.write('cont_mat\n')
    for cm in con_mat:
        f.write(cm+'\n')
    f.close()

if __name__ == "__main__":
    # Simple test for atoms 
    # creat input 
    write_pyfodmc_atoms(sys='Kr')
    # Fortran call 
    fodmc.fodmc.points_on_sphere_metropolis_spin_centers()

    # Simple test for molecule
    sys = 'SO2.xyz'
    con_mat = ['(1-2)-(2-2)','(1-3)-(2-2)','1-(1-1)','2-(2-2)','3-(2-2)']
    write_pyfodmc_molecules(sys=sys,con_mat=con_mat)
    # Fortran call 
    fodmc.fodmc.points_on_sphere_metropolis_spin_centers()
