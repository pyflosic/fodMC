from fodMC.pyfodmc import pyfodmc
from fodMC.pyfodmc.mol2fodmc import mol2fodmc
import time

molecule = 'C2H4'
mol2fodmc(molecule+'.mol')
t1 = time.time()
pyfodmc.get_guess(molecule)
#pyfodmc.get_guess()
t2 = time.time()
print('fodMC FOD generation for {}: {} s'.format(molecule,t2-t1))
