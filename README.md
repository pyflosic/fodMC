# fodMC - Fermi-orbital descriptor Monte-Carlo 
[![license](https://img.shields.io/badge/license-APACHE2-green)](https://www.apache.org/licenses/LICENSE-2.0)
[![language](https://img.shields.io/badge/language-Fortran90-red)](https://www.fortran90.org/)
[![language](https://img.shields.io/badge/language-Python3-blue)](https://www.python.org/)
[![version](https://img.shields.io/badge/version-1.0.14-lightgrey)](https://github.com/pyflosic/fodMC/blob/master/README.md)
[![doi](https://img.shields.io/badge/DOI-10.5281/zenodo.3922473-blue)](https://zenodo.org/record/3922473#.Xvor5aYpCCg)

Main developer:  

*  Kai Trepte (KT, kai.trepte1987@gmail.com)    
*  Sebastian Schwalbe (SS, theonov13@gmail.com)  

Sidekicks:  
  
* Alex Johnson (AJ, johns1ai@cmich.edu)   
* Jakob Kraus (JaK, jakob.kraus@physik.tu-freiberg.de)   

## Description
fodMC is a generator for Fermi-orbital descriptor (FOD) positions      
to be used in the Fermi-Löwdin orbital self-interaction correction (FLO-SIC) method.               
There is a publication, explaining the underlying idea of this program; please see    


### [Interpretation and automatic generation of Fermi-orbital descriptors](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26062)
```
S. Schwalbe, K. Trepte, et al.      
Journal of Computational Chemistry, vol. 40, pp. 2843-2857, 2019    
```

## Installation (using pip) 

```bash 
$ pip3 install fodMC
```


## Installation (local)

```bash 
$ git clone https://github.com/pyflosic/fodMC
$ cd fodMC
$ pip3 install -e .
```

The Python module is called fodmc. 

There are examples to make you familiar with the execution.    
You can either work at the Fortran level (see examples/fortran) or    
use the Python interface, see examples/python.   


# ATTENTION
Initial FODs for transition metals and larger atoms may currently not be reliable.       
This is due to the fact that the spherical symmetry of the core FODs (with is strictly enforced)     
might not represent a good guess for such systems. Intensive research is needed to determine the     
correct structural motifs for such systems.     
