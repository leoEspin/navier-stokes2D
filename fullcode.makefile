#the main code:
mcode=solute.half
#extra code, the solver:
excode1=pSolver
excode2=exactSoln.half
excode3=nsSolver
excode4=adi_solute.half
#libraries, extra flags:
comp=gfortran -Wall -O3 -ffast-math -fopenmp#-msse2 -mfpmath=sse  #-std=f95
libs=-lfftw   
exec=-o oxc.full
addr=-L.
main: $(mcode).o $(excode2).o $(excode3).o $(excode4).o
	$(comp) $(exec) $(mcode).o $(excode1).o $(excode2).o $(excode3).o $(excode4).o $(addr) $(libs)
$(mcode).o: $(mcode).f90  $(excode3).o $(excode4).o
	$(comp) -c $(mcode).f90
$(excode1).o: $(excode1).f90
	$(comp) -c $(excode1).f90
$(excode2).o: $(excode2).f90
	$(comp) -c $(excode2).f90
$(excode3).o: $(excode3).f90 $(excode1).o $(excode2).o
	$(comp) -c $(excode3).f90
$(excode4).o: $(excode4).f90 $(excode3).o $(excode2).o
	$(comp) -c $(excode4).f90      
