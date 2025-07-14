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



# two approaches to computing the physical time at iom
# kk = (theta_at_ioc - d - theta0[i])/p
# k = round(kk)
#
# Now we compute the physical time t at which the buffer occupancy is measured.
# We could instead use the approach: physical_time_at_iom = invert(theta[i], theta[i].y[end] - d)
# which would be equally good.
# we know theta0[i] + k*p = theta_at_ioc - d
#
# physical_time_at_iom = invert_with_guess(theta[i], theta0[i] + k*p, guess2)
@inline function get_occupancy(physical_time, theta, theta_at_dst, ugn, edge, latency, gear, guess1)
    theta_at_src_delayed_by_latency = evaluate_with_guess(theta[edge.src], physical_time - latency, guess1)
    return ugn + floor(gear*theta_at_src_delayed_by_latency) - floor(gear*theta_at_dst)
end

# measurement
function set_measurement_values(i, status)
    # This means in particular that sampling always starts
    # off roughly in phase across all nodes. Perhaps we should have
    # an option where we can test different initial phasing.
    #
    # That is, the first measurement at node i always
    # happens at physical time t = 0, and so the first sample
    # must happen at physical time t =  d/omega[i]_{-1}
    #
    #
    theta_at_ioc = status.theta[i].y[end]
    theta_at_iom = theta_at_ioc - status.model.d
    physical_time_at_iom = invert_with_guess(status.theta[i], theta_at_iom, status.guess2)
    for (j,e) in enumerate(status.model.graph.incoming_edges[i])
        mm = get_occupancy(physical_time_at_iom, status.theta, theta_at_iom, status.model.ugns[e],
            status.model.graph.edges[e], status.model.latency[e], status.model.gears[e], status.guess1)
        status.measurement[i].occupancies[j] = mm
        log_occupancy(status.model.slog, e, physical_time_at_iom, mm)
    end
    status.measurement[i].theta_at_ioc = theta_at_ioc
    status.measurement[i].physical_time_at_iom = physical_time_at_iom
end


function make_initial_state(i, epoch, freq_m2, freq_m1, theta0, d, slog)
    sample_times = Float64[epoch,                             0]
    theta_values = Float64[theta0 + epoch * freq_m2,          theta0]
    if d > 0
        push!(sample_times, d/freq_m1)
        push!(theta_values, theta0 + d)
    end
    theta_i = PiecewiseLinear(sample_times, theta_values)
    return theta_i
end

function make_guesses(errors, n, latency, d)
    wavg = sum(initial_error.(errors)) / n
    latency_min = minimum(latency)
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

function find_node_whose_freq_changes_next(theta)
    # The following is equivalent to
    #   minimizing_physical_time, minimizing_nodeid  = findmin([theta[j].x[end] for j=1:n])
    # but is a faster version
    n = length(theta)
    minimizing_nodeid = 1
    minimizing_physical_time = theta[1].x[end]
    for i=2:n
        t = theta[i].x[end]
        if t < minimizing_physical_time
            minimizing_physical_time = t
            minimizing_nodeid = i
        end
    end
    return Event(minimizing_nodeid, minimizing_physical_time)
end

function all_nodes_finished(theta, t)
    for th in theta
        if th.x[end] < t
            return false
        end
    end
    println("C1: all nodes finished")
    return true
end

function make_initial_status(model::Model)
    n = model.graph.n
    incoming = model.graph.incoming_edges
    measurement = [Measurement(zeros(length(incoming[i]))) for i = 1:n]
    guess1, guess2 = make_guesses(model.errors, n, model.latency, model.d)
    initial_state(i) = make_initial_state(i, model.epoch, model.wm2[i], model.wm1[i], model.theta0[i], model.d, model.slog)
    theta = [initial_state(i) for i=1:n]
    status = Status(; model, measurement, theta, guess1, guess2)
    return status
end

function tick(status, i)
    set_measurement_values(i, status)
    c::Float64 = next(status.model.controller[i], status.measurement[i])
    return c
end

# event.time is the physical time of next frequency change
# event.node is the node id
function update(event, status)
    i = event.node
    t = event.time
    c::Float64 = tick(status, i)
    # If the error is constant, then ds = p / (c + e.const_error)
    # This is the physical time between the current control action and the next one.
    # That is, how much physical time must pass, starting from time t, until the clock increases by p?
    physical_duration_until_next_freq_change = local_to_realtime(status.model.errors[i], status.model.p, c,
        t, status.model.wmin)
    check_freq_is_positive(status.model.checker, c, status.model.errors, i, t)
    extend_dx_dy(status.theta[i], physical_duration_until_next_freq_change, status.model.p)
    done = sim_is_done(status.model.stopper, status.model.errors, i, status.model.p, c, t)
    return done
end

# We cannot have theta0 = 0. Because frames are popped from the elastic
# buffer and sent when theta is 0,1,2,3,4.
# The controller makes measurements and updates frequency when
# theta is theta0 + k*poll_period
# If theta0 is an integer then the controller is evaluating theta
# at times when it an integer, and so it may evaluate it as k - epsilon
# or k + epsilon. These cases are different, since in one case
# a frame has been sent, and in the other it has not.




function runcallisto(model::Model)
    n = model.graph.n
    check_initial(model.p, initial_error.(model.errors), n, model.tmax)
    status = make_initial_status(model)
    done = false
    while (!all_nodes_finished(status.theta, status.model.tmax)) && !done
        event = find_node_whose_freq_changes_next(status.theta)
        done = update(event, status)
    end
    return model.slog, status.theta
end

function callisto(c)
    print("Running callisto1: ")
    compute_parameters(c)
    @time simlog, theta = runcallisto(c.model)
    return (simlog = simlog, theta = theta,)
end










end
