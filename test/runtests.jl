
module RunTests

include("testsim.jl")
using .TestSim

include("testopts.jl")
using .TestOpts

using Test

function main()
    @testset "TestSim" begin
        @test TestSim.main1()
        @test TestSim.main2()
        @test TestSim.main3()
    end
    @testset "TestOpts" begin
        @test TestOpts.main1()
        @test TestOpts.main2()
        @test TestOpts.main3()
    end
end


end

using .RunTests
RunTests.main()



