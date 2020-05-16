from ase.atoms import Atoms

def write_fodMC_system(ase_atoms,cm,options={'unit' : 'angstrom','fix1s' : 'fix1s'}):
    ''' Writes fodMC system file '''
    f = open('system','w') 
    sym = ase_atoms.get_chemical_symbols() 
    natoms = len(sym) 
    pos = ase_atoms.get_positions() 
    f.write('{}\n'.format(natoms))
    str_o = ''
    for o in list(options.keys()):
        str_o += options[o] + ' '
    f.write(str_o[:-1]+'\n')
    for a in range(natoms):
        f.write('{} {} {} {} \n'.format(sym[a],pos[a][0],pos[a][1],pos[a][2]))
    f.write('con_mat\n') 
    f.write(cm)
    f.write('\n\n')
    f.close() 

def read_atoms_bond_mol(f_name):
    ''' Adjusted ase MDF mol (chemical table format) reader''' 
    # Input:  mol file 
    # Output: ase_atoms, connectivity matrix (cm) 
    # Notes:  16.05.2020 -- currently only supports single, double, trible bonds 
    #                       one need to descided how to input 5-3 etc bonds 

    fileobj = open(f_name,'r') 
    lines = fileobj.readlines()
    L1 = lines[3]
    
    # The V2000 dialect uses a fixed field length of 3, which means there
    # won't be space between the numbers if there are 100+ atoms, and
    # the format doesn't support 1000+ atoms at all.
    if L1.rstrip().endswith('V2000'):
        natoms = int(L1[:3].strip())
        nbonds = int(L1[3:6].strip())
    else:
        natoms = int(L1.split()[0])
        nbonds = int(L1.split()[1])
    positions = []
    symbols = []
    for line in lines[4:4 + natoms]:
        x, y, z, symbol = line.split()[:4]
        symbols.append(symbol)
        positions.append([float(x), float(y), float(z)])
    # Bonding types 
    BOtype = {'1' : '(1-1)','2' : '(2-2)','3' : '(3-3)'}
    # Connectivity matrix 
    cm = '' 
    for l in range(4+natoms,4+natoms+nbonds):
        line = lines[l]
        A, B, BO = line.split()[:3]
        cm +='({}-{})-{}'.format(A,B,BOtype[BO])+' '
    return Atoms(symbols=symbols, positions=positions),cm 

def mol2fodmc(mol):
    ''' Converts mol2 format to fodMC system '''
    ase_atoms, cm = read_atoms_bond_mol(mol)
    write_fodMC_system(ase_atoms=ase_atoms,cm=cm)

