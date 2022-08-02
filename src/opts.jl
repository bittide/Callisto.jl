
module Opts

export CalOpts

using Random

import ..Topology: Graph
import ..SimCore: Error


"""
   Link structure

"""
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

"""
    CalOpts structure

"""
mutable struct CalOpts
    graph
    links
    tmax
    kp
    ki
    poll_period
    control_delay
    theta0
    errors
    wm1
    wm2
    wmin
    epoch
    controller_init
    controller_next
end



# assume f(a) < 0 and f(b) > 0
function bisection(f, a, b, tol = 1e-5)
    while true
        c = (a + b)/2
        if abs(a-b) < tol
            return c
        end
        y = f(c)
        if y == 0
            return c
        end
        if y < 0
            a = c
        else
            b = c
        end
    end
end


"""
    a tuplmake_ugn(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst) 

Compute the value for the UGN of a particular link, given

 - `beta0`: the initial occupancy of the elastic buffer at the destination
 - `gear`: the gear ratio of the link, that is the number of frames sent per source clock tick
 - `latency`: the time taken for a frame to traverse the link
 - `theta0_at_src`: the value of the clock phase at time 0 at the source end of the link
 - `wm2_at_src`: the clock frequency at the source end of the link for times < 0
 - `theta0_at_dst`: the value of the clock phase at time 0 at the destination end of the link

"""
function make_ugn(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst)
    return beta0 - floor(gear*(theta0_at_src - latency * wm2_at_src)) + floor(gear*theta0_at_dst)
end


function make_frequencies(seed, num_nodes)
    rng = Random.Xoshiro(seed)
    f = 1 .+  rand(rng, 1:10^5, num_nodes) *  1e-9
    return f
end


    
"""
    CalOpts(;)

Return a CalOpts structure. Keyword arguments are:

   graph
   latency
  
"""
function CalOpts(;topology = ("mesh", 3, 2), graph = nothing,
                 kp = 2e-8, ki = 1e-15, tmax = 1e9,
                 latency=5000, seed=1, theta0=0.1,
                 poll_period = 1e5,
                 control_delay = 1000,
                 bidirectional = true,
                 base_freq = 1,
                 wmin = 0.1,
                 wm1 = nothing, wm2 = nothing,
                 controller_init = nothing, controller_next = nothing,
                 errors=nothing, gears = nothing, beta0 = 50)

    if !isnothing(graph)
        g = graph
    else
        g = Graph(topology; bidirectional)
    end
    
    if isa(latency, Number)
        latency = latency*ones(g.m)
    end
    
    epoch = -10*maximum(latency)

    # errors is a function which takes
    # node i, local time interval p, correction c, and wall-clock time s
    # and returns a wall-clock time interval ds
    #
    # so that int_s^{s+ds} (c + \omega_i(t)) dt = p
    #
    # If omega_i is constant, then
    #
    #   ds*(c + omega_i) = p
    #
    # that is, ds = p / (c + omega_i)
    #
    if isnothing(errors)
        errors_const = make_frequencies(seed, g.n)
        errors =  [Error(errors_const[i]) for i=1:g.n]
    end


    if isa(wm1, Number)
        wm1 = wm1 * ones(g.n)
    elseif isnothing(wm1)
        wm1 = [ errors[i](0) for i=1:g.n]
    end


    if isa(wm2, Number)
        wm2 = wm2 * ones(g.n)
    elseif isnothing(wm2)
        wm2 = [ errors[i](0) for i=1:g.n]
    end



    # instead of beta[dst,src]
    #
    # use beta[edgeid] where src,dst = g.edges[edgeid]
    #
    if isa(beta0, Number)
        beta0 = beta0 * ones(g.m)
    end
   
    
    measurement_offsets = beta0 

    if isa(theta0, Number)
        theta0 = theta0*ones(g.n)
    end

    if isnothing(gears)
        gears = ones(g.m)
    end

    function make_link(e)
        return Link(e, g.edges[e].src, g.edges[e].dst,
                    make_ugn(beta0[e], gears[e], latency[e], theta0[g.edges[e].src],
                        wm2[g.edges[e].src], theta0[g.edges[e].dst]),
                    latency[e],
                    gears[e],
                    beta0[e],
                    measurement_offsets[e])
    end

    links = [make_link(e) for e = 1:g.m]
    
    if isnothing(controller_init)
        controller_init = () -> 0.0
    end

    function cnext(i, xi, measurement)
        if length(measurement) == 0
            # can have no measurements if
            # there are no incoming edges at a node
            # and hence no elastic buffers
            return 0,0
        end
        r = sum(a[2] for a in measurement)
        next_state =  xi + poll_period*r/base_freq
        correction =  kp * r +  ki * next_state
        return next_state, correction
    end
    
    if isnothing(controller_next)
        controller_next = cnext
    end


    c = CalOpts(g,
                links,
                tmax,
                kp,
                ki,
                poll_period,
                control_delay,
                theta0,
                errors,
                wm1,
                wm2,
                wmin,
                epoch,
                controller_init,
                controller_next
                )

    return c
end

# 
# function nscale(x)
#     s = repr(x)
#     if 'e' in s
#         s = "\\texttt{$s}"
#     end
#     return s
# end
# 
# function xprintln(c::CalOpts)
#     topology = c.graph.topology
#     if topology[1] == "full" && topology[2] == 2
#         s = "2 nodes"
#     elseif topology[1] == "full" && topology[2] == 3
#         s = "triangle"
#     else
#         s = topology[1] * "(" * join(string.(topology[2:end]),", ") *")"
#     end
#     println("\\item topology: $s")
#     println("\\item \$k_P =  $(nscale(c.kp))\$,  \$k_I = $(nscale(c.ki))\$")
# 
#     n = c.graph.n
#     latency = c.links[1].latency
#     # warning: assumes uniform latency
#     println("\\item \$\\codestyle{latency} = $(latency)\$")
# 
#     println("\\item \$\\codestyle{poll\\_period} = $(c.poll_period)\$")
# 
#     #s = "(" * join(repr.(c.errors), ", ") * ")"
#     #println("\\item \$\\codestyle{uncorrected_frequency} = $s\$")
#     println("\\item \$\\codestyle{control\\_delay} = $(c.control_delay)\$")
#     println("\\item \$\\codestyle{theta0} = $(c.theta0)\$")
# end
# 



end

