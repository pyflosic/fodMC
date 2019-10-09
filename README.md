![GitHub Logo](images/logo_fodMC.png)

# fodMC - Fermi-orbital descriptor Monte-Carlo 

Main developer:  

*  Kai Trepte (KT, trept1k@cmich.edu)    

Sidekicks:  

* Alex Johnson (AJ, johns1ai@cmich.edu)   
* Jakob Kraus (JaK, jakob.kraus@student.tu-freiberg.de)   
* Sebastian Schwalbe (SS, theonov13@gmail.com)    

Coding language: FORTRAN   

## Description
   
Guess generator for FOD positions to be used in the Fermi-Löwdin orbital self-interaction correction (FLO-SIC) method.           
There is a manual, explaining the underlying idea of this program and showing some examples for its usage.

## Installation 
The fodMC is written in FORTRAN. Make sure you have a fortran compiler like gfortran to compile the code.  
To compile the code, go to the *src* directory and type   

```compile
make
```
```alternative debricated script to compile
        bash compile.sh
```

There are now two executables. 'fodMC' is the original code, while 'fodMC_motif' uses structural motifs 
for the generation of atomic guesses and core FODs. Thus, in this version these FODs are not generated 
by the distribution of points on a sphere (beta test version).

# ATTENTION
Initial FODs for transition metals and larger atoms may currently not be reliable.
This is due to the fact that the spherical symmetry of the core FODs (with is strictly enforced)
might not represent a good guess for such systems. Intensive research is needed to determine the 
correct structural motifs for such systems.


## Running the code 

Go to the folders *examples*. 

You might want to change the relative directory paths in the run.sh to your absolute src-directory (see manual as well).                                         
```Usage bash
        bash run.sh
```
```should work with all shells
./run.sh
```

to run the code for a given input in the file *system*.
If you want to use e.g. 'fodMC_motif' instead of 'fodMC', 
you need to change the executable name in the run.sh file.

## Versions
* The original fodMC code is called *fodMC_alex.f90*. 
* The next major update was included in *fodMC_motif.f90*, where structural motifs are used for core FODs (instead of using MC to get them).

## Manual/Tutorial

There is a manual in the 'doc' folder, called fodMC.pdf. 
Furthermore, there are several tutorial videos about the 
usage and handling of the fodMC at youtube, 
at the channel 'The extended Physiker Clan'.
