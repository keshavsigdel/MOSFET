FC=gfortran -c
LD=gfortran
SRC=linspace.f90 mysum.f90 fermi.f90 main.f90
OBJ=linspace.o mysum.o fermi.o main.o

fig111:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o fig111.x
	rm -rf $(OBJ)
clean:
	rm -rf fig111.x $(OBJ)
