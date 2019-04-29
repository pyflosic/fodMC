# PyFODMC - A simple python backend for the fodMC 
Coding language: python  
Author: Seabstian Schwalbe  
Licence: The same as the fodMC, because it is part of it. 

* compile.sh:   compiling fodmc as libary (changes in fodmc_alex.f90 needed)
                
* pyfodmc.py:  use fodmc.so as python module 
* src fodmc:    the src is adjusted module header/subroutine 

## Installation 
You need the fodmc and f2py compiler. The f2py compiler is included in the numpy package.
```bash 
bash compile.sh
```

## How to use 
Simply execute the pyfodmc.py script in the bash. There are some examples 
included. (This should be extented in the future.) 

```bash 
python pyfodmc.py
```

To clean the directory, use the following command 
```bash
bash clean.sh 
```
