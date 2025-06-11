
module Sim

export callisto
import ..SimCore: SimCore, Link, Model, Cal1Opts, compute_parameters, log_occupancy
import Piecewise: PiecewiseLinear, evaluate_with_guess, extend_dx_dy, invert, invert_with_guess
import StateSystems: next
import Topology: inport_indexed_from_edge_indexed

using BittideTools
using JuliaTools
using Random


##############################################################################



##############################################################################

@inline get_phase(theta::PiecewiseLinear, t::Float64, guess::Int) = evaluate_with_guess(theta, t, guess::Int)
#
@inline function  betax((;ugn, src, dst, latency, gear), t, theta, theta_at_iom, guess1)
    phase1 = get_phase(theta[src], t - latency, guess1)
    return ugn + floor(gear*phase1) - floor(gear*theta_at_iom)
end

@inline function get_measurement(link, physical_time_at_iom, theta, theta_at_iom, guess1)
    return betax(link, physical_time_at_iom, theta, theta_at_iom, guess1)
end
##############################################################################

# measurement
function Gx(i, theta, s, theta0, p, links, incoming, d, slog, measurement, guess1, guess2)
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

    indegree = length(incoming[i])
    occ = measurement.occupancies
    for (j,e) in enumerate(incoming[i])
        link = links[e]
        mm = get_measurement(link, physical_time_at_iom, theta, theta_at_iom, guess1)
        occ[j] = mm # (mm - link.offset)/link.gear
        log_occupancy(slog, e, physical_time_at_iom, mm)
    end
    measurement.theta_at_iofc = theta_at_iofc
    measurement.physical_time_at_iom = physical_time_at_iom
end


function initial_statex(i, epoch, freq_m2, freq_m1, theta0, d, slog)
    sample_times = Float64[epoch,                             0,                 d/freq_m1]
    theta_values = Float64[theta0 + epoch * freq_m2,          theta0,            theta0 + d]
    theta_i = PiecewiseLinear(sample_times, theta_values)
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
    wavg = sum(initial_error.(errors)) / num_nodes
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
                  stopper,
                  controller, checker)

    check(d, p, initial_error.(errors), num_nodes, tmax)
    measurement = [Measurement(zeros(length(incoming[i]))) for i = 1:num_nodes]
    guess1, guess2 = make_guesses(errors, num_nodes, links, d)
    initial_state(i) = initial_statex(i, epoch, freq_m2[i], freq_m1[i], theta0[i], d, slog)
    theta = [initial_state(i) for i=1:num_nodes]
    s = 0.0
    done = false
    while s < tmax && !done
        #s, i  = findmin([theta[j].x[end] for j=1:num_nodes])
        s, i  = find_soonest_node(theta, num_nodes)

        Gx(i, theta, s, theta0, p,
           links, incoming, d,
           slog, measurement[i],
           guess1, guess2)
        #
        # constants:        i, indegree
        # immutable params: theta_at_iofc, physical_time_at_iom
        # mutable refs:     occ_measurements
        #
        #measurements = (; occ_measurements, theta_at_iofc, physical_time_at_iom)
        c::Float64 = next(controller[i], measurement[i])

        # If the error is constant, then ds = p / (c + e.const_error)
        # This is the physical time between the current control action
        # and the next one.
        ds = local_to_realtime(errors[i], p, c, s, wmin)
        if !check_freq_is_positive(checker, c, errors, i, s)
            break
        end
        extend_dx_dy(theta[i], ds, p)
        done = sim_is_done(stopper, errors, i, p, c, s)
    end
    return slog, theta
end

function runcallisto(model::Model)
    return callisto(model.slog, model.tmax, model.epoch, model.num_nodes, model.links,
        model.incoming, model.p, model.d, model.errors, model.theta0, model.wm2,
        model.wm1, model.wmin, model.stopper, model.controller, model.checker)
end




function callisto(c)
    print("Running callisto: ")

    compute_parameters(c)

    # callisto cares about only the incoming edges to each node
    @time simlog, theta = runcallisto(c.model)

    return (simlog = simlog, theta = theta,)
end










end
