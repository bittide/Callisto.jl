
module Callisto

include("simcore.jl")
using .SimCore

# depends on SimCore and Topology 
include("opts.jl")
using .Opts

# depends on SimCore, LogData, and Piecewise
include("sim.jl")
using .Sim

# depends on SimCore, LogData, and Piecewise
include("post.jl")
using .Post

# depends on post
include("fastpostproc.jl")
using .FastPostProc


function reexport(m)
    for a in names(m)
        eval(Expr(:export, a))
    end
end

reexport(SimCore)
reexport(Opts)
reexport(Sim)
reexport(Post)
reexport(FastPostProc)


end

