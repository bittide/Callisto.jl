
module Opts

export CalOpts

using Random

import ..Topology: Graph
import ..SimCore: Error, initial_error, SimCore


"""
   Link structure

"""
mutable struct Link
    id::Int64
    src::Int64
    dst::Int64
    ugn::Int64
    latency::Float64
    gear::Float64
    beta0::Int64
    offset::Float64
end


mutable struct Report
    avgfreq
    niters
    expected_log_length
    expected_time
    latencies
    ugns
    beta0
    initial_frames_on_link
end

function Base.show(io::IO, r::Report)
    f(x) = round(x; sigdigits = 2)
    g(x, y) = println(io, rpad(x, 30), "= ", y)

    g("expected log length", r.expected_log_length)
    g("expected run time", f(r.expected_time))
    g("latencies", f.(r.latencies))
    g("ugns", r.ugns)
    g("initial buffer occupancies", r.beta0)
    g("initial link occupancies", r.initial_frames_on_link)
end

function Report(graph, errors, poll_period, tmax, links, theta0, wm2)
    n = graph.n
    m = graph.m
    # expected simulation time
    avgfreq = sum(initial_error.(errors))/n
    niters = Int.(round.(avgfreq * n * tmax / poll_period))
    expected_log_length = niters + 3*n
    expected_time = expected_log_length * 5e-6
    latencies = [l.latency for l in links]
    ugns = [Int(l.ugn) for l in links]
    beta0 = [Int(l.beta0) for l in links]

    initial_frames_on_link = [Int(get_initial_frames_on_link(links[e].gear,
                                                         links[e].latency,
                                                         theta0[graph.edges[e].src],
                                                         wm2[graph.edges[e].src]))
                              for e=1:m]
                                                         
    
    r = Report(avgfreq,
               niters,
               expected_log_length,
               expected_time,
               latencies,
               ugns,
               beta0,
               initial_frames_on_link)
    return r
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
    beta_maker
    base_freq
    special
    stopper
    report
    safe
    enable_logging
    controller
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


function get_initial_frames_on_link(gear, latency, theta0_at_src,  wm2_at_src)
    return  floor(gear*theta0_at_src) - floor(gear*(theta0_at_src - latency * wm2_at_src)) 
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
                 beta_maker = nothing,
                 wm1 = nothing, wm2 = nothing,
                 controller = nothing,
                 errors=nothing, gears = nothing, beta0 = 50,
                 safe = true,
                 enable_logging = true,
                 special = 0)

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
        errors = make_frequencies(seed, g.n)
        #errors =  [Error(errors_const[i]) for i=1:g.n]
    end


    if isa(wm1, Number)
        wm1 = wm1 * ones(g.n)
    elseif isnothing(wm1)
        wm1 = [ initial_error(errors[i]) for i=1:g.n]
    end


    if isa(wm2, Number)
        wm2 = wm2 * ones(g.n)
    elseif isnothing(wm2)
        wm2 = [ initial_error(errors[i]) for i=1:g.n]
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

    function cnext_sum(i, xi, measurement)
        r = measurement
        next_state =  xi + poll_period*r/base_freq
        correction =  kp * r +  ki * next_state
        return next_state, correction
    end


    stopper = nothing

    report = Report(g, errors, poll_period, tmax, links, theta0, wm2)

    #if isnothing(controller)
    #    controller =  [SimCore.PIController(kp, ki, 0.0, poll_period, base_freq) for i=1:g.n]
    #end
                
                

    
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
                beta_maker,
                base_freq,
                special,
                stopper,
                report,
                safe,
                enable_logging,
                controller
                )

    return c
end



end

