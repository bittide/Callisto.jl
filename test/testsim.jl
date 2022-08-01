

module TestSim

##############################################################################
module Piecewise

struct PiecewiseLinear 
    x
    y
end

interpolate(x1,y1,x2,y2,d) = ((x2-d)*y1 - (x1-d)*y2)/(x2-x1)

function evaluatepl(x, y, t)
    # this function returns the largest i such that x[i] <= t
    i = searchsortedlast(x, t)
    if x[i] == t
        return y[i]
    end
    @assert i != length(x)
    return interpolate(x[i], y[i], x[i+1], y[i+1], t)
end

invert(p::PiecewiseLinear, w) = evaluatepl(p.y, p.x, w)
(p::PiecewiseLinear)(t) = evaluatepl(p.x, p.y, t)

end


##############################################################################
module LogData

Log(a,b) = return
llog(a,b,c) =  return
llognext(a) =  return

end


##############################################################################
module SimCore
beta((;ugn, src, dst, latency, gear), t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))


struct ConstError
    const_error::Float64
end

local_to_realtime(e::ConstError, p, c, s, wmin) = p / (c + e.const_error)
(e::ConstError)(t) = e.const_error

end


##############################################################################
struct Link
    id
    src
    dst
    ugn
    latency
    gear
    beta0
    offset
end


##############################################################################

# what we are actually testing
include("../src/sim.jl")
using .Sim

function makelink(e, src, dst, beta0, theta0, latency, wm2, gear)
    ugn = beta0 - floor(gear*(theta0[src] - latency * wm2[src])) + floor(gear*theta0[dst])
    return Link(e, src, dst, ugn, latency, gear, beta0, beta0)
end

function get_incoming_edges(edges, n)
    incoming_edges = [Int64[] for i=1:n]
    for (i, (src, dst)) in enumerate(edges)
        push!(incoming_edges[dst], i)
    end
    return incoming_edges
end

function offsetnormsquared(p::Piecewise.PiecewiseLinear, offset)
    s = 0.0
    for i=1:length(p.x)-1
        dt = (p.x[i+1] - p.x[i])
        dtheta = p.y[i+1] - p.y[i]
        freq = dtheta/dt
        s += (freq - offset)^2 * dt
    end
    return s
end

function testcallisto(tmax,
                      epoch,
                      num_nodes,
                      edges,
                      poll_period,
                      d,
                      kp,
                      ki,
                      freqs)


    theta0 = [0.1, 0.1, 0.1]
    errors = [SimCore.ConstError(a) for a in freqs]
    wm2 = copy(freqs)
    wm1 = copy(freqs)
    num_edges = length(edges)
    beta0 = fill(50, num_edges)
    latency = 5000
    links = [makelink(e, edges[e][1], edges[e][2], beta0[e], theta0, latency, wm2, 1) for e = 1:length(edges)]
    incoming = get_incoming_edges(edges, num_nodes)
    wmin = 0.1
    
    controller_init = () -> 0.0
    function controller_next(i, xi, ml)
        s = sum(a[2] for a in ml)
        next_state =  xi + poll_period*s
        correction =  kp * s +  ki * next_state
        return next_state, correction
    end

    slog, theta = callisto(tmax,
                           epoch,
                           num_nodes,
                           links,
                           incoming,
                           poll_period,
                           d,
                           errors,
                           theta0, 
                           wm2,
                           wm1,
                           wmin,
                           controller_init,
                           controller_next)

    return theta
end


function main()
    kp = 2e-8
    ki = 1e-15
    tmax = 200e6
    epoch = -5e4
    num_nodes = 3
    alpha = 1e-4
    f0 = 1
    freqs = [f0 + alpha, f0-alpha, f0]
    edges = ( (1,2), (1,3), (2,3), (2,1), (3,1), (3,2))
    poll_period = 100000
    d = 1000
    theta =  testcallisto(tmax,
                          epoch,
                          num_nodes,
                          edges,
                          poll_period,
                          d,
                          kp,
                          ki,                          
                          freqs)

    @assert all([a.x[end] > tmax for a in theta]) "Error: Callisto simulation did not run for required interval."
    ns = sum(offsetnormsquared(q, f0) for q in theta)

    # resistance distance test
    res = 2/3
    expected_ns = alpha^2*res/(2*kp)
    @assert  0.9 <= ns/expected_ns <= 1.1 "Error: Callisto simulation frequency not within predicted bounds."

    return true
end
    


end
