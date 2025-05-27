
module TestPiecewise



using Piecewise
using Test


##############################################################################


function main1()
    x = [1,2,4,7,8]
    y = [ 1,3,1,3 ]
    p = PiecewiseConstant(x,y)
    @assert tuples(p)[1] == (1,1)
    @assert tuples(p)[2] == (2,1)
    @assert tuples(p)[3] == (2,3)
    @assert tuples(p)[4] == (4,3)
    @assert tuples(p)[5] == (4,1)
    @assert tuples(p)[6] == (7,1)
    @assert tuples(p)[7] == (7,3)
    @assert tuples(p)[8] == (8,3)
    @assert length(tuples(p)) == 8
    return true
end

function main2()
    x = [1,2,4,7]
    y = [1,3,1,3]
    p = PiecewiseLinear(x,y)
    @assert tuples(p)[1] == (1,1)
    @assert tuples(p)[2] == (2,3)
    @assert tuples(p)[3] == (4,1)
    @assert tuples(p)[4] == (7,3)
    @assert length(tuples(p)) == 4
    return true
end

function main3()
    x = [1,2,4,7]
    y = [1,3,1,3]
    p = Samples(x,y)
    @assert tuples(p)[1] == (1,1)
    @assert tuples(p)[2] == (2,3)
    @assert tuples(p)[3] == (4,1)
    @assert tuples(p)[4] == (7,3)
    @assert length(tuples(p)) == 4
    return true
end

function main4()
    x = [1,2,4,12]
    y = [1,3,1,3]
    p = PiecewiseLinear(x,y)
    @assert p(1) == 1
    @assert p(2) == 3
    @assert p(3) == 2
    @assert p(6) == 1.5
    return true
end

function main99()
    return true
end


function main()
    @testset "Callisto piecewise" begin
        @test main1()
        @test main2()
        @test main3()
        @test main4()

        @test main99()
    end
    return true
end



end

