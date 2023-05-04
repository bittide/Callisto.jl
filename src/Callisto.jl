
module Callisto

include("piecewise.jl")
using .Piecewise

include("logdata.jl")
using .LogData

include("topology.jl")
using .Topology

include("simcore.jl")
using .SimCore

include("statedynamics.jl")
using .StateDynamics

# depends on SimCore and Topology 
include("opts.jl")
using .Opts

# depends on SimCore, LogData, and Piecewise
include("sim.jl")
using .Sim

# depends on SimCore, LogData, and Piecewise
include("post.jl")
using .Post


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
export beta, gamma
#------------------------------------------------------------------------------
# exports from opts.jl

export CalOpts

#------------------------------------------------------------------------------
# exports from sim.jl

export callisto

#------------------------------------------------------------------------------
# exports from post.jl

export parse_callisto_log, parse_callisto_logx, focused_callisto_info, get_freq
export getlatency, getlambda, getbeta, getomega, getdata





end

