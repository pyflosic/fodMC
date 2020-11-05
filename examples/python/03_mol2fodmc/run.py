from fodMC.pyfodmc import pyfodmc
from fodMC.pyfodmc.mol2fodmc import mol2fodmc
import glob 
import time 
p = './structures/'
mol = glob.glob(p+'*.mol')
for m in mol:
    mol2fodmc(m)
    tmp_name = str(m).split('/')[-1].split('.')[0]
    t1 = time.time()

    output_mode = 'PyFLOSIC'
    output_name = m.split('/')[-1].replace('.mol','.xyz')
    pyfodmc.get_guess(output_mode,output_name)

    t2 = time.time() 
    print('fodMC FOD generation for {}: {} s'.format(tmp_name+'.mol',t2-t1))


