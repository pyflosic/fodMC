# fodMC      

Fermi-orbital descriptor Monte-Carlo 
Main developer: Kai Trepte      
Coding language: FORTRAN   
                                         
   
Guess generator for FOD positions to be used in the 
Fermi-Löwdin self-interaction correction (FLO-SIC) method           
#  

There is a manual, explaining the underlying idea of this program and showing some examples for its usage.

The fodMC is written in FORTRAN. Make sure you have a fortran compiler like gfortran to compile the code.
To compile the code, go to the 'src' directory and type   

        bash compile.sh


Go to the folders 'examples'. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!! Change the directory in the run.sh to your src-directory !!!   
!!! See manual as well                                       !!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
Then, type 

        bash run.sh

to run the code for a given input in the file 'system'.
