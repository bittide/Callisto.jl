
module TestStateDynamics

include("../src/statedynamics.jl")

using .StateDynamics: StateSystem, PIStateSystem, next, StoppingStateSystem, SequentialStateSystem
using Test

apply(K::StateSystem, inputs) = [next(K, u) for u in inputs]


# Verify that the integrator integrates
function main1()
    inputs = [1.0, 2, 3, 7]
    K = PIStateSystem(0, 2)
    outputs = apply(K, inputs)
    @assert outputs == 2 * cumsum(inputs)
    return true
end

# Verify that proportional controllers just scale inputs
function main2()
    inputs = [1.0, 2, 3, 7]
    K = PIStateSystem(3, 0)
    outputs = apply(K, inputs)
    @assert outputs == 3 * inputs
    return true
end

# Verify that PI controllers are P + I 
function main3()
    inputs = [1.0, 2, 3, 7]
    Kp = PIStateSystem(3, 0)
    Ki = PIStateSystem(0, 2)
    Kpi = PIStateSystem(3, 2)
    outp = apply(Kp, inputs)
    outi = apply(Ki, inputs)
    outpi = apply(Kpi, inputs)
    @assert outpi == outp + outi
    return true
end

# Verify that a stopping system stops
function main4()
    koff = 4
    Kpi = PIStateSystem(3, 2)
    K = StoppingStateSystem(koff, Kpi)
    inputs1 = [1.0, 2, 3, 7, 8, 4, 2, 1]
    # inputs2[t] = 0 for t > koff 
    # The input values at those times do not affect the output
    # because by then the system has stopped.
    inputs2 = inputs1 .* [(i <= koff) for i = 1:length(inputs1)]
    outputs1 = apply(K, inputs1)

    # We recreate system to reset the state to zero
    Kpi = PIStateSystem(3, 2)
    K = StoppingStateSystem(koff, Kpi)
    outputs2 = apply(K, inputs2)
    @assert outputs1 == outputs2
    return true
end

# Verify that a switching system switches
function main5()
    ksw = 4
    K1 = PIStateSystem(3, 2)
    K2 = PIStateSystem(1, 3)
    swf(lastoutput) = return
    K = SequentialStateSystem([K1, K2], [ksw], [swf])
    inputs1 = [1.0, 2, 3, 7, 8, 4, 2, 1]
    outputs1 = apply(K, inputs1)

    # reset by recreating
    K1 = PIStateSystem(3, 2)
    K2 = PIStateSystem(1, 3)
    # Run the two systems seprately, one after the other
    outputsfirst = apply(K1, inputs1[1:ksw])
    outputssecond = apply(K2, inputs1[ksw+1:end])

    @assert outputs1 == vcat(outputsfirst, outputssecond)

    return true
end

    
    



    
function main()
    @testset "Callisto statedynamics" begin
        @test main1()
        @test main2()
        @test main3()
        @test main4()
        @test main5()
    end
    return true
end



end
