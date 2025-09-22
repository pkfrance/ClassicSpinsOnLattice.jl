module ClassicSpinsOnLattice

using StaticArrays: SVector, SMatrix, deleteat, insert, pushfirst

export Link, @Link, SpinProblem, add_link!, unit_vector, du_dp

include("problem.jl")
include("spin.jl")

end
