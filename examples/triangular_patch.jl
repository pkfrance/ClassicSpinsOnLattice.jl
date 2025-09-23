using ClassicSpinsOnLattice

J=-1.0
sp=SpinProblem{2}(2) # Create Hamiltonian for 2D spins on a 2D lattice
add_link!(sp, @Link((1,1), [1, 0], J))
add_link!(sp, @Link((1,1), [0, 1], J))
add_link!(sp, @Link((1,1), [1, 1], J))

