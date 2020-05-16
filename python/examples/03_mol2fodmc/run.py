from pyfodmc import pyfodmc
from pyfodmc.mol2fodmc import mol2fodmc

p = './structures/'
mol = [p+'291.mol',p+'6085.mol',p+'6086.mol',p+'236.mol']
for m in mol:
    mol2fodmc(m)
    tmp_name = str(m).split('/')[-1].split('.')[0]
    pyfodmc.get_guess(name=tmp_name)


