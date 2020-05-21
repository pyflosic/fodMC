![GitHub Logo](fodMC/doc/images/logo_fodMC.png)

# fodMC - Fermi-orbital descriptor Monte-Carlo 

Main developer:  

*  Kai Trepte (KT, trept1k@cmich.edu)    

Sidekicks:  

* Sebastian Schwalbe (SS, theonov13@gmail.com)    
* Alex Johnson (AJ, johns1ai@cmich.edu)   
* Jakob Kraus (JaK, jakob.kraus@physik.tu-freiberg.de)   

Coding language: FORTRAN, Python 

## Description
The fodMC is a generator for Fermi-orbital descriptor (FOD) positions to be used in the Fermi-Löwdin orbital self-interaction correction (FLO-SIC) method.           
There is a publication, explaining the underlying idea of this program; please see

'Interpretation and automatic generation of Fermi-orbital descriptors', S. Schwalbe, K. Trepte, et al., Journal of Computational Chemistry, vol. 40, pp. 2843-2857, 2019


## Installation 
The fodMC is written in FORTRAN. Make sure you have a fortran compiler like gfortran to compile the code.  
To compile the code, go to the *lib* directory and type   

```compile
make
```

There are now two executables. 'fodMC' is the original code, while 'fodMC_motif' uses structural motifs 
for the generation of atomic guesses and core FODs. Thus, in this version these FODs are not generated 
by the distribution of points on a sphere.

# ATTENTION
Initial FODs for transition metals and larger atoms may currently not be reliable.
This is due to the fact that the spherical symmetry of the core FODs (with is strictly enforced)
might not represent a good guess for such systems. Intensive research is needed to determine the 
correct structural motifs for such systems.


## Running the code (in the shell, no GUI) 
Go to the folders *examples/fortran*. 

You might want to change the relative directory paths in the run.sh to your absolute src-directory.  
```should work with all shells
./run.sh
```

to run the code for a given input in the file *system*.
If you want to use e.g. 'fodMC' instead of 'fodMC_motif', 
you need to change the executable name in the run.sh file.

## Running the code (Python-based GUI)
Go to the folder *python_GUI*.
In the *examples* folder, one can find how to use the graphical user interface (GUI) of the fodMC. 

Thanks goes to Sebastian Schwalbe for developing this GUI.

## Running the code (Python)
In the folder *pyfodMC*, there is a Python-Overlay of the fodMC. 

Thanks goes to Sebastian Schwalbe to develop this.

## Tutorials

There are several tutorial videos about the 
usage and handling of the fodMC at youtube, 
at the channel 'The extended Physiker Clan'.

In addition, there is a publication summarizing essential findings regarding FODs, see 
'Interpretation and automatic generation of Fermi-orbital descriptors', S. Schwalbe, K. Trepte, et al., Journal of Computational Chemistry, vol. 40, pp. 2843-2857, 2019
