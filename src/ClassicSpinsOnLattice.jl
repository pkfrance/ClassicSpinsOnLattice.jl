module ClassicSpinsOnLattice

using StaticArrays: SVector, deleteat, insert, pushfirst

export Link, @Link, SpinProblem, add_link!, dir

include("problem.jl")
include("spin.jl")

end
