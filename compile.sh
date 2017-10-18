# Clearing
binDir="bin/" ; comando="cd $binDir"; eval $comando
comando="rm *.o"; eval $comando
comando="rm *.mod"; eval $comando

# Compiling modules
comando="gfortran -c ../src/termo.F90"; eval $comando
comando="gfortran -c ../src/meshStructure.F90"; eval $comando
comando="gfortran -c ../src/scalarStructure.F90"; eval $comando
comando="gfortran -c ../src/solver.F90"; eval $comando
comando="gfortran -c ../src/io.F90"; eval $comando
comando="gfortran -c ../src/mInputReader.F90"; eval $comando
comando="gfortran -c ../src/shapeFunctions.F90"; eval $comando
comando="gfortran -c ../src/setup.F90"; eval $comando
comando="gfortran -c ../src/utilities.F90"; eval $comando
comando="gfortran -c ../src/scalar.F90"; eval $comando
comando="gfortran -c ../src/postProc.F90"; eval $comando
comando="gfortran -c ../src/processor.F90"; eval $comando
# Linking and compiling driver
comando="gfortran ../src/driver.F90 termo.o meshStructure.o scalarStructure.o solver.o io.o mInputReader.o setup.o shapeFunctions.o utilities.o  scalar.o postProc.o processor.o -o ../finel.exe"; eval $comando
