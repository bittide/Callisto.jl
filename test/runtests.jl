
module TestSet


include("testsim.jl")
using .TestSim

include("testpiecewise.jl")
using .TestPiecewise

using Test

@testset "Callisto tests" begin
    @test TestSim.main()
    @test TestPiecewise.main()
end
    

end


