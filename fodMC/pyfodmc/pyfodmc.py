# pyfodmc - python interface for fodmc code to pyflosic 
# Author:       Sebastian Schwalbe (SS)  
# Changelog:    11.02.2019 
#               SS update to the new fodmc.f90 
#               SS add new write_pyfodmc_molecules function 
#               15.05.2020 
#               re-included the python interface 
import fodmc 
from fodMC.pyfodmc.database import write_database
from ase.io import read,write  
import os

def get_database(input_file='xx_database_xx'):
    f = open(os.path.dirname(os.path.abspath(__file__))+'/'+input_file,'r') 
    ll = f.readlines()
    f.close() 
    p = os.getcwd() 
    o = open(p+'/'+input_file,'w')
    for l in range(len(ll)):
        o.write(ll[l]) 
    o.close() 

def clean_files(): 
    files = ['CLUSTER','FRMORB','system','xx_database_xx']
    for f in files: 
        try:
            if os.path.exists(f):
                os.remove(f)
        except: 'Nothing' 

def rename_xyz(name):
    os.rename('Nuc_FOD.xyz',name+'.xyz')


def get_guess(name='fodMC'):
    # cp the database to the calculation folder 
    #get_database(input_file='xx_database_xx')
    write_database()
    # magic to capture that output:
    # from http://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python
    #      http://websrv.cs.umt.edu/isis/index.php/F2py_example
    output_file = name + '.out' 
    if os.path.exists(output_file):
        os.remove(output_file)
    # open outputfile
    outfile = os.open(output_file, os.O_RDWR|os.O_CREAT)
    # save the current file descriptor
    save = os.dup(1)
    # put outfile on 1
    os.dup2(outfile, 1)
    # end magic
    # FORTAN call 
    fodmc.fodmc_mod.get_guess()
    # restore the standard output file descriptor
    os.dup2(save, 1)
    # close the output file
    os.close(outfile)
    f = open(output_file,'r')
    output = f.read()
    f.close()
    # rm files which are not needed 
    clean_files()
    # rename output xyz
    rename_xyz(name)

def write_pyfodmc_atoms(sys):
    #
    # write fodmc input for atoms 
    #
    # create system file 
    f = open('system','w')
    f.write('1 %s\n' % sys)
    f.write('angstrom fix1s\n')
    f.write('%s 0.0 0.0 0.0\n\n' %sys)
    f.close()

def write_pyfodmc_molecules(sys,con_mat):
    
    ase_atoms = read(sys)
    sym = ase_atoms.get_chemical_symbols()
    pos = ase_atoms.get_positions() 
    natoms = len(ase_atoms.get_chemical_symbols())
    
    # create system file 
    f = open('system','w')
    f.write('%i %s\n' % (natoms,sys))
    f.write('angstrom fix1s\n')
    for p in range(len(pos)):
        f.write('%s %0.5f %0.5f %0.5f\n' %(sym[p],pos[p][0],pos[p][1],pos[p][2]))
    f.write('cont_mat\n')
    for cm in con_mat:
        f.write(cm+'\n')
    f.close()

if __name__ == "__main__":
    def make_atom():
        # Simple test for atoms 
        # creat input 
        write_pyfodmc_atoms(sys='Kr')
        # Fortran call 
        fodmc.fodmc.get_guess()
    
    def make_molecule():
        # Simple test for molecule
        sys = 'SO2.xyz'
        con_mat = ['(1-2)-(2-2)','(1-3)-(2-2)\n']   
        write_pyfodmc_molecules(sys=sys,con_mat=con_mat)
        # Fortran call 
        fodmc.fodmc.get_guess()
    make_atom()
