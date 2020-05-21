.SUFIXES: .f90

# Compiler
FC = gfortran

# Flags
# Flags for debugging 
#FFLAGS = -Wall -O0
# Flags for production
FFLAGS = -O3
P1=fodMC
P2=fodMC_motif

# Default target

all: fodMC fodMC_motif

fodMC:  $(P1)_alex.f90
	$(FC) $(P1)_alex.f90 -o $(P1) $(FFLAGS)

fodMC_motif: $(P2).f90
	$(FC) $(P2).f90 -o $(P2) $(FFLAGS)

clean:
	rm -f fodMC fodMC_motif 
