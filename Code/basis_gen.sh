#This shell script will generate the basis
#states needed to run the hof_solver_2.f90 program
#It will delete any file called out.dat though!

#ARGS: N, Q
gfortran -O3 -o bstate bstate.f90
FIL="./out.dat"
if [ -e "$FIL" ]
then
	ls
	rm out.dat
fi
./bstate $1 $2
mv out.dat inff1.dat
