from pyfodmc import pyfodmc
from pyfodmc.mol2fodmc import mol2fodmc
import glob 
import time 
p = './structures/'
mol = glob.glob(p+'*.mol')
for m in mol:
    mol2fodmc(m)
    tmp_name = str(m).split('/')[-1].split('.')[0]
    t1 = time.time() 
    pyfodmc.get_guess(name=tmp_name)
    t = time.time() - t1 
    print('{}: {} s'.format(tmp_name+'.mol',t))


