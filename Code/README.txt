INTERACTIONS:

(See NB:1, pg 170 for derivation of interactions term)

This program bundle's purpose is to create and solve
a Hamiltonian in the k space for the tight binding
model, with a magnetic field perpendicular to the
plane, for an arbitary number of particles.
Now, with interactions as well.
it is basically the same as hof_solver_2.f90,
but now has on the diagonal an additional
term that is proportional to the N
Overview:
solve_BH.sh is the over-arching script.
First if checks to make sure that the directory
that you are working in is 'clear'. Ie, that you
don't have certain keyword files in your directory.
Then, it will create a temporary directory and
move into it.

This is done so that you can ./solve_BH.sh more
than one time, and not have an issue with
cross contamination of data files/inputs/temperary
files created during the process.

Then, it enters into the basis state creator/checker.
It checks to see if you have a directory in your
working directory called BasisStates, and if you do,
then it additionally checks to see if you have the
right basisstate file.
If you do, then it will move the correct basisstate file
into the temporary directory, and call the fortran ho_solver_2.f90

If you don't, then the program will create a directory called BasisStates, and then call the fortran
function bstate.f90, which will then create the basis
states using an algorithm found in a paper by Carl de Boor (1999)
"Computational aspects of multivariate polynomial interpolation: Indexing the coefficients"

It will then save this file in the folder 
BasisStates.

NOTE********************
This was originally done as, naively, I though that
we would be reaching N=100, Q=8, where the calculation
of basis states is non-trivial, and if you had already 
calculated the basis states for a N,Q, you could
save a lot of computational time.

However, as we are working in the N<8 regime, it is
trivially calculated.
END NOTE*********************

Basis Encoding Schema:
The basis states are in k space.
They are encoded as an Q long list. N is
the total number of particles.
Each element represents a k state, and so, it 
can take on N number of values.
The trick is that N must be conserved, so each row
must add up to N.
Example:
N = 2, Q = 3

           0           0           2
           0           1           1
           0           2           0
           1           0           1
           1           1           0
           2           0           0

So the first row shows all two particles
being created in the k(3) state.
The second row shows one particle in the k(2) state
and one particle in the k(3) state.

After the basis state issue is sorted out, then the 
main program is called: hof_solver_2.f90.

HOF_SOLVER_INT.f90.

Now with INTERACTIONS!

The main change in this program from all the others,
is that it directly uses the basis states to 
calculate the appropriate prefactors, allowed
transitions, etc.
So for each element in the H matrix H(B(i),B(j)), where
B(i) and B(j) denote the basis states), we need
to calculate the appropriate value.

It was noted that for a transition between
basis states to occur, there can only be SINGLE
difference seperating them. 
So, to see if an element in the hamiltonian will
be non-zero, we can take two rows and minus them.
This is the purpose of the 'transtest' array
in the code. It stores the result of this minus.

If a transition between the two basis states is allowed, 
then we can predict that the resulting array stored in
transtest will contain:
one -1
one 1
the rest of the elements will be zero.

Example:
To test for a transition between B(1) and B(2):
B(1) - B(2) = [0 , -1, 1]

Further more, one can think of these -1, 1's as
creation and anhilation operators that link
these two basis states.

B(1) <-> B(5) 
B(1) - B(5) = [-1, -1, 2]

This test fails, as it should.

So, if the transition is allowable.

However, we have further information about the
Hamiltonian that we can use.

We know that the only elements in the Hamiltonians
equation that allow transitions between states are the
exp(-i*ky) and exp(i*ky) terms, and they ONLY allow 
transitions between neighbor states.
More specifically, this means if we represent
the creation and anhilation operator as above-
as a list containing 1 for creation, and -1 for
anhilation - then the only allowed transitions
are the ones where the 1 and the -1 are neighbors
on the list. Ie [0,0,0-1,1,0] would be allowed
but [1,0,0,0,-1] would not be.

This knowledge allows us to further prune those
states that will be zero in the H matrix.

The final piece of knowledge used is the order in
which the creation and anhilation operators are in.

When I derived the Harper equations using the
creation/anhilation operators, the two exp terms looked
like:

1) exp(-i*ky)*bstar(k-1)*b(k)
2) exp(i*ky)*bstar(k+1)*b(k)


This is the final piece of information that allows
us to exactly determine the hamiltonian.

If transtest = [...,-1,1,...] then we
know that the term that goes there is exp(-i*ky)
mutiplied by the appropriate prefactor.

If transtest = [...1,-1...] then we know that
the term that goes there is exp(i*ky) multiplied
by the appropriate prefactor.

Once the hamiltonian is determined, it is a simple
matter to feed it into zheev and solve for the
eigenvalues.
These are then put into a folder called Res (for
Results). The output of the hof_solver_2.f90 is
placed in a folder called Logs.

The script may give you some warnings about:
"There is no such file" but you can ignore them.


