
function dir(params::SVector{1, Float64})
    return SVector{2, Float64}((cos(params[1]), sin(params[1])))
end

function dir(params::SVector{D, Float64}) where D
    tail=sin(params[1])*dir(deleteat(params, 1))
    pushfirst(tail, cos(params[1]))
end