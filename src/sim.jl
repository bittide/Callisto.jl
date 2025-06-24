# Copyright 2023 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

module Sim

using ..SimCore
using BittideTools
using JuliaTools
using Piecewise
using Random
using StateSystems
using Topology

export callisto

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

    sample_times = Float64[epoch,                             0]
    theta_values = Float64[theta0 + epoch * freq_m2,          theta0]
    if d > 0
        push!(sample_times, d/freq_m1)
        push!(theta_values, theta0 + d)
    end

    theta_i = PiecewiseLinear(sample_times, theta_values)
    return theta_i
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


# We cannot have theta0 = 0. Because frames are popped from the elastic
# buffer and sent when theta is 0,1,2,3,4.
# The controller makes measurements and updates frequency when
# theta is theta0 + k*poll_period
# If theta0 is an integer then the controller is evaluating theta
# at times when it an integer, and so it may evaluate it as k - epsilon
# or k + epsilon. These cases are different, since in one case
# a frame has been sent, and in the other it has not.
#
function callisto(slog, tmax, epoch, num_nodes, links, incoming, p, d, errors,
                  theta0, freq_m2, freq_m1, wmin,
                  stopper,
                  controller, checker)

    check_initial(p, initial_error.(errors), num_nodes, tmax)
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
    print("Running callisto1: ")
    compute_parameters(c)
    @time simlog, theta = runcallisto(c.model)
    return (simlog = simlog, theta = theta,)
end










end
