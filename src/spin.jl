

"""
    unit_vector(p::SVector{D, Float64}) where D

Returns a unit vector of dimension `D+1` parameterized by `D` angles. 
For `D=1` and `D=2`, the parameterization corresponds to the classical 
polar and spherical coordinates, respectively. The input `p` is a vector 
of angles.

# Type Parameters
- `D`: Dimension of the input angle vector.

# Arguments
- `p::SVector{D, Float64}`: Vector of angles used to parameterize the unit vector.

# Returns
- `SVector{D+1, Float64}`: A unit vector in `D+1` dimensions.
"""
function unit_vector(p::SVector{D, Float64}) where D
    tail=sin(p[1])*unit_vector(deleteat(p, 1))
    pushfirst(tail, cos(p[1]))
end

"""
    unit_vector(p::SVector{1, Float64})
Specialized version of `unit_vector` for 2D vectors.
"""
function unit_vector(p::SVector{1, Float64})
    return SVector{2, Float64}((cos(p[1]), sin(p[1])))
end

"""
    du_dp(p::SVector{D, Float64}) where D
Computes the Jacobian matrix of the `unit_vector` function with respect to the input angles `p`.
# Type Parameters
- `D`: Dimension of the input angle vector.
# Arguments
- `p::SVector{D, Float64}`: Vector of angles used to parameterize the unit vector.
# Returns
- `SMatrix{D+1, D, Float64}`: Jacobian matrix of the `unit_vector` function.
"""
function du_dp(p::SVector{D, Float64}) where D
    b=sin(p[1])*du_dp(deleteat(p, 1))
    r=vcat(zeros(Float64, D-1), b)
    v=cos(p[1])*unit_vector(deleteat(p, 1))
    c=pushfirst(v, -sin(p[1]))
    hcat(c, r)
end

"""
    du_dp(p::SVector{1, Float64})
Specialized version of `du_dp` for 2D vectors.
"""
function du_dp(p::SVector{1, Float64}) 
    return SMatrix{2, 1, Float64}((-sin(p[1]), cos(p[1])))
end