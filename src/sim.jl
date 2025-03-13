


module Sim

#using DataStructures
#using InteractiveUtils


export callisto


struct Checker end


sim_is_done(stopper::Nothing, errors, i, p, c, s) = false

import ..Piecewise:  PiecewiseLinear, invert, evaluate_with_guess, invert_with_guess
import ..LogData:    Log, llog, llognext, LogData
import ..SimCore:    beta, local_to_realtime, errorat, SimCore, initial_error, init, getlog, step, Controller

# IOM = instant of measurement
# IOFC = instant of frequency change

# also needed in opts.jl
# This is useful for very large simulations
# where we give up flexibility and logging for speed

##############################################################################

function clog(slog::Log, sk, nid, freq)
    llog(slog, 1, 1)
    llog(slog, 3, sk)
    llog(slog, 4, nid)
    llog(slog, 6, freq)
    llognext(slog)
end

function clog(slog::Log, k, sk, nid, freq)
    llog(slog, 1, 1)
    llog(slog, 2, k)
    llog(slog, 3, sk)
    llog(slog, 4, nid)
    llog(slog, 6, freq)
    llognext(slog)
end

clog(slog::Nothing, sk, nid, freq) = return
clog(slog::Nothing, k, sk, nid, freq) = return
LogData.llog(x::Nothing, a, b) = return

@inline log_controller(slog::Nothing, controller, occ_measurements, i, indegree) = return 

function log_controller(slog::Log, controller, occ_measurements, i, indegree)
    llog(slog, 7, sum_or_missing(a for a in occ_measurements[1:indegree]))
    llog(slog, 8, getlog(controller))
end

log_w(slog::Nothing, errors, i, p, c, s, wmin, theta) = return

function log_w(slog::Log, errors, i, p, c, s, wmin, theta)
    w = c + errorat(errors[i], s)    # only used for log and termination
    clog(slog, theta[i].x[end], i, w)
end

@inline loop_check(ch::Nothing, c, errors, i, s) = return

function loop_check(ch::Checker, c, errors, i, s)
    if c + errorat(errors[i], s) <= 0
        println("Frequency at node $i has become negative: omega = ", c + errors[i](s))
        println("time s = ", s)
    end
end
##############################################################################

@inline get_phase(theta::PiecewiseLinear, t::Float64, guess::Int) = evaluate_with_guess(theta, t, guess::Int)
#
@inline function  betax((;ugn, src, dst, latency, gear), t, theta, theta_at_iom, guess1)
    phase1 = get_phase(theta[src], t - latency, guess1)
    return ugn + floor(gear*phase1) - floor(gear*theta_at_iom)
end


@inline function get_measurement(b::Nothing, link, physical_time_at_iom, theta, theta_at_iom, guess1)
    return (betax(link, physical_time_at_iom, theta, theta_at_iom, guess1) - link.offset)/link.gear
end
##############################################################################

# measurement
function Gx(i, theta, s, theta0, p, links, incoming, d, beta_maker, slog, occ_measurements, guess1, guess2)
    #
    # This means in particular that sampling always starts
    # off roughly in phase across all nodes. Perhaps we should have
    # an option where we can test different initial phasing.
    #
    # That is, the first measurement at node i always
    # happens at physical time t = 0, and so the first sample
    # must happen at physical time t =  d/omega[i]_{-1}
    #
    #
    theta_at_iofc = theta[i].y[end]
    theta_at_iom = theta_at_iofc - d
    k = round((theta_at_iofc - d - theta0[i])/p)
    # Now we compute the physical time t at which the buffer occupancy is measured.
    # We could instead use the approach: physical_time_at_iom = invert(theta[i], theta[i].y[end] - d)
    # which would be equally good.
    # we know theta0[i] + k*p = theta_at_iofc - d
    physical_time_at_iom = invert_with_guess(theta[i], theta0[i] + k*p, guess2)
    llog(slog, 2, k)
    llog(slog, 5, physical_time_at_iom)
    indegree = length(incoming[i])
    for (j,e) in enumerate(incoming[i])
        occ_measurements[j] = get_measurement(beta_maker, links[e], physical_time_at_iom,
                                   theta, theta_at_iom, guess1)
    end
    return theta_at_iofc, physical_time_at_iom, indegree
end

function sum_or_missing(x)
    if length(x) == 0
        return -100
    end
    return sum(x)
end

# This should be preallocated, since we know tmax
@inline function extend(p::PiecewiseLinear, dx, dy)
    push!(p.x, p.x[end] + dx)
    push!(p.y, p.y[end] + dy)
end

@inline function Fx(i, theta, ds, p)
    # p/(c+e) is the time interval between s_k and s_{k+1}
    extend(theta[i], ds, p)
end

function initial_statex(i, epoch, freq_m2, freq_m1, theta0, d, slog)
    sample_times = Float64[epoch,                             0,                 d/freq_m1]
    theta_values = Float64[theta0 + epoch * freq_m2,          theta0,            theta0 + d]
    theta_i = PiecewiseLinear(sample_times, theta_values)
    clog(slog, -2, epoch, i, freq_m2)
    clog(slog, -1, 0, i, freq_m1)
    return theta_i
end

function check(d, p, w, n, tmax)
    # I don't see any theoretical or practical reason
    # for the restriction d < p, so this has been removed.
    avgfreq = sum(w)/n
    if avgfreq*tmax/p < 10
        println("Warning; fewer than 10 steps per node in simulation")
    end
    relmadf = sum(abs.(w .- avgfreq))/avgfreq
    if relmadf < 1e-9
        println("Warning: frequencies are very close together. Expect strange behavior.")
    end
end

function make_guesses(errors, num_nodes, links, d)
    wavg = sum(SimCore.initial_error.(errors)) / num_nodes
    latency_min = minimum([l.latency for l in links])
    # The two values guess1 and guess2 are approximately how many ticks ago
    # we need to evaluate theta in order to compute buffer occupancy.
    # If we are wrong, it's better to be slighly too small.
    # But the code will work if these are wrong; they simply
    # allow speedup heuristics.
    guess1 = Int(round(d + wavg * latency_min - 1))
    guess2 = Int(round(d - 1))
    if guess1 < 0
        guess1 = 0
    end
    if guess2 < 0
        guess2 = 0
    end
    return guess1, guess2
end

function find_soonest_node(theta, n)
    minind = 1
    minval = theta[1].x[end]
    for i=2:n
        s = theta[i].x[end]
        if s < minval
            minval = s
            minind = i
        end
    end
    return minval, minind
end

function callisto(slog, tmax, epoch, num_nodes, links, incoming, p, d, errors,
                  theta0, freq_m2, freq_m1, wmin, 
                  special, stopper, maxindegree,
                  controller, checker, beta_maker)

    check(d, p, initial_error.(errors), num_nodes, tmax)
    occ_measurements = zeros(maxindegree)
    guess1, guess2 = make_guesses(errors, num_nodes, links, d)
    initial_state(i) = initial_statex(i, epoch, freq_m2[i], freq_m1[i], theta0[i], d, slog)
    theta = [initial_state(i) for i=1:num_nodes]

    s = 0.0
    done = false
    while s < tmax && !done
        #s, i  = findmin([theta[j].x[end] for j=1:num_nodes])
        s, i  = find_soonest_node(theta, num_nodes)
        theta_at_iofc, physical_time_at_iom, indegree = Gx(i, theta, s, theta0, p,
                                                           links, incoming, d, beta_maker,
                                                           slog, occ_measurements,
                                                           guess1, guess2)
        log_controller(slog, controller[i], occ_measurements, i, indegree)
        c::Float64 = step(controller[i], i, theta_at_iofc,
                          physical_time_at_iom, indegree, occ_measurements)
        # If the error is constant, then ds = p / (c + e.const_error)
        # This is the physical time between the current control action
        # and the next one.
        ds = local_to_realtime(errors[i], p, c, s, wmin)
        log_w(slog, errors, i, p, c, s, wmin, theta)
        loop_check(checker, c, errors, i, s)
        Fx(i, theta, ds, p)
        done = sim_is_done(stopper, errors, i, p, c, s)
    end
    return slog, theta
end



function callisto(c)
    print("Running callisto: ")

    if c.enable_logging
        header = ["ltype", "k", "sk", "nid", "tk", "freq", "meas", "xi"]
        slog = Log(8,header)
    else
        slog = nothing
    end
    if c.safe
        checker = Checker()
    else
        checker = nothing
    end


    maxindegree = maximum([length(c.graph.incoming_edges[i]) for i = 1:c.graph.n])


    
    # callisto cares about only the incoming edges to each node
    @time simlog, theta = callisto(slog, c.tmax, c.epoch, c.graph.n, c.links, c.graph.incoming_edges, c.poll_period, c.control_delay, c.errors, c.theta0, c.wm2, c.wm1, c.wmin, c.special, c.stopper, maxindegree, c.controller, checker, c.beta_maker)

    return (simlog = simlog, theta = theta,)
end










end
