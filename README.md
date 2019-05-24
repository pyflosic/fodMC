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

        bash compile.sh

There are now two executables. 'fodMC' is the original code, while 'fodMC_test' uses structural motifs 
for the generation of atomic guesses and core FODs. Thus, in this version these FODs are not generated 
by the distribution of points on a sphere (beta test version).

## Running the code 

Go to the folders *examples*. 

You might want to change the relative directory paths in the run.sh to your absolute src-directory (see manual as well).                                         

Then, type 

        bash run.sh

to run the code for a given input in the file *system*.
If you want to use 'fodMC_test', you need to change the executable name in the run.sh file.

## Manual/Tutorial

There is a manual in the 'doc' folder, called fodMC.pdf. 
Furthermore, there are several tutorial videos about the 
usage and handling of the fodMC at youtube, 
at the channel 'The extended Physiker Clan'.
