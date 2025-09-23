using Test
using ClassicSpinsOnLattice
using LinearAlgebra: norm, dot
using StaticArrays: SVector

@testset "unit_vector norm test" begin
    for d in 1:4
        for _ in 1:20
            p = 2π .* rand(SVector{d, Float64})
            v = unit_vector(p)
            @test isapprox(norm(v), 1.0; atol=1e-15)
        end
    end
end

@testset "jacobian test" begin
    for d in 1:4
        for _ in 1:10
            p = 2π .* rand(SVector{d, Float64})
            v = unit_vector(p)
            J = du_dp(p)
            # Test shape
            @test size(J) == (d+1, d)
            # Test that J columns are orthogonal to v
            for i in 1:size(J, 2)
                @test isapprox(dot(v, J[:, i]), 0.0; atol=1e-14)
            end
            for __ in 1:10
                dp=1e-8*(rand(SVector{d, Float64}).-0.5)
                dv1=unit_vector(p+dp)-v
                dv2=J*dp
                @test isapprox(norm(dv1-dv2), 0.0; atol=1e-15)
            end
        end
    end
end
