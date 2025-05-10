
module Callisto

using Piecewise
using LogData
using Topology
using StateSystems


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



#------------------------------------------------------------------------------
# re-exports from piecewise

export tuples
export PiecewiseConstant, PiecewiseLinear, Series, Samples
export integer_crossings, evaluate, differentiate, 
    floor, *, +, delay, integrate, after, before, IntegersBetween, stairs,
    xy, limits, crossing, xmax, invert, discontinuities, square, l2norm,
    steady_state, remove_steady_state, integer_crossings_in_interval,
    definite_integral


#------------------------------------------------------------------------------
# exports from Topology
export resistance

#------------------------------------------------------------------------------
# exports from simcore.jl

export Error
export beta, gamma, errorat, initial_error, controllerfunctions
#------------------------------------------------------------------------------
# exports from opts.jl

export CalOpts, set_controller

#------------------------------------------------------------------------------
# exports from sim.jl

export callisto

#------------------------------------------------------------------------------
# exports from post.jl

export parse_callisto_logx, focused_callisto_info


#------------------------------------------------------------------------------
# exports from fastpost

export  FastPost, fast_post_process, quick_post, callistox


end

