#!/bin/sh

# Author: S. Schwalbe 
# Date: 27.09.2019 
# you might need to change the paths according to your needs
fodmc="$(dirname "$(pwd)")"
echo $fodmc
cp $fodmc/src/xx_database_xx . 
#$fodmc/src/fodMC
$fodmc/src/fodMC_motif
