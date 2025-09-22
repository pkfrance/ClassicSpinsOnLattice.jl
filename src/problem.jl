
"""
    Link(i1, i2, offset, J)

Represents a link between two sites with spins in a lattice.

# Arguments
- `i1::Int`: Index of the first site (must be strictly positive).
- `i2::Int`: Index of the second site (must be strictly positive).
- `offset`: Lattice translation vector linking the unit cells of sites `i1` and `i2`. Should be zero if both sites belong to the same unit cell.
- `J`: Interaction term. `J > 0` indicates anti-ferromagnetic interaction, while `J < 0` indicates ferromagnetic interaction.
"""
struct Link{N}
    i1::Int
    i2::Int
    offset::SVector{N,Int}
    J::Float64

    # Internal constructor with validation
    function Link{N}(i1::Int, i2::Int, offset::SVector{N,Int}, J::Float64) where N
        if i1 <= 0 || i2 <= 0
            throw(ArgumentError("Site indices i1 and i2 must be strictly positive. Got i1=$i1, i2=$i2."))
        end
        new(i1, i2, offset, J)
    end
end


"""
    @Link((i1, i2), offset, J)

Convenience macro to construct a `Link` object. 
- `(i1, i2)`: Tuple of site indices (both strictly positive).
- `offset`: Vector of integers representing the lattice translation.
- `J`: Interaction term.

Automatically infers the lattice dimension `N` from the length of `offset`.
"""
macro Link(indices, offset, J)
    # Extract indices from the tuple
    i1, i2 = esc.(indices.args)
    # Convert offset to SVector and get dimension
    offset_expr = esc(offset)
    J_expr = esc(J)

    return quote
        local _N = length($offset_expr)
        local _offset_vec = SVector{_N}($offset_expr)

        Link{_N}($i1, $i2, _offset_vec, $J_expr)
    end
end

"""
    struct SpinProblem{N}

Represents the Hamiltonian problem for a classical spin system on a lattice of dimension `N`, 
where each spin has `dim` components. The lattice contains `n_sites` per unit cell. 
The `links` vector specifies the topology of connections between the sites, defining 
the interaction structure of the lattice.

# Arguments
- `dim::Int`: Number of components for each classical spin.
- `N::Int`: Dimension of the lattice.
- `n_sites::Int`: Number of sites per unit cell.
- `links::Vector{...}`: Describes the connectivity between sites.

# Type parameters
- `N::Int`: Dimension of the lattice.

"""
mutable struct SpinProblem{N}
    links::Vector{Link{N}}
    n_sites::Int
    dim::Int
end
"""
    SpinProblem{N}(dim::Int)

Constructs an empty `SpinProblem` for a lattice of dimension `N` and spins with `dim` components.
The `links` vector is initialized empty and `n_sites` is set to zero.

# Arguments
- `dim::Int`: Number of components for each classical spin.

# Returns
A `SpinProblem{N}` instance with no links and zero sites.
"""
function SpinProblem{N}(dim::Int) where N
    SpinProblem{N}(Link{N}[], 0, dim)
end

"""
    add_link!(sp::SpinProblem{N}, l::Link{N})

Adds a link (interaction term) to the Hamiltonian of the spin problem `sp`.
If either site index in `l` exceeds the current number of sites per unit cell,
`n_sites` is updated accordingly.

# Arguments
- `sp::SpinProblem{N}`: The spin problem to modify.
- `l::Link{N}`: The link to add.

# Returns
The modified `SpinProblem{N}` (in-place).
"""
function add_link!(sp::SpinProblem{N}, l::Link{N})::SpinProblem{N} where N
    push!(sp.links, l)
    sp.n_sites = max(sp.n_sites, l.i1, l.i2)
    return sp
end
