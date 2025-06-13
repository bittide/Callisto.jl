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

module SimCore

using Topology
using JuliaTools
using BittideTools
using LogData

export Link, CalOpts, Cal1Opts, compute_parameters, log_init, log_occupancy

const log_fields = (; tk=1, mocc=2)
log_init(m) =  [Log(log_fields) for i=1:m]
@inline function log_occupancy(slog::Vector{Log}, eid, t, mocc)
    LogData.log(slog[eid], log_fields.tk, t)
    LogData.log(slog[eid], log_fields.mocc, mocc)
    start_new_record(slog[eid])
end
@inline log_occupancy(slog::Nothing, eid, t, mocc) = return


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

Base.@kwdef mutable struct Model
    slog
    tmax
    epoch
    num_nodes
    links
    incoming
    p
    d
    errors
    theta0
    wm2
    wm1
    wmin
    stopper
    controller
    checker
end

Base.@kwdef mutable struct ComputedParams
    output_systems
    offset
    local_offsets
    latency
end


"""
    CalOpts structure
Holds configuration options for a Callisto simulation.

# Fields
- `graph`:                      The network topology. Defaults to a 3-node bidirectional triangle.
- `latency`:                    Link latencies. Can be a scalar (applied to all links) or a vector. Defaults to 20.
- `offset`:                     Initial buffer offsets/target occupancies. Can be a scalar or a vector. Defaults to 200.
- `tmax`:                       Maximum simulation time. Defaults to 1.5e5.
- `kp`:                         Proportional gain for the PI controller. Defaults to 8e-6.
- `ki`:                         Integral gain for the PI controller. Defaults to 0.0.
- `poll_period`:                Time between controller updates for each node. Defaults to 2000.
- `theta0`:                     Initial phase for each node. Can be a scalar or a vector. Defaults to 0.0.
- `errors`:                     Clock frequency errors for each node.
                                Can be `:random` (generates random errors), a vector of numerical error values,
                                or a vector of `SimCore.Error` objects. Defaults to `:random`.
- `base_freq`:                  Base frequency used by the controller for scaling the integral term. Defaults to 1.0.
- `stopper`:                    A custom stopping condition for the simulation. Defaults to `nothing` (no custom stopper).
- `check_freq_is_positive`:     If `true`, checks if node frequencies become negative during simulation. Defaults to `true`.
- `enable_logging`:             If `true`, enables logging of simulation data (e.g., buffer occupancies). Defaults to `true`.
- `controller`:                 Specifies the controller type for each node.
                                Can be `:pi` (default PI controller) or a vector of
                                custom controller objects. Defaults to `:pi`.
- `v`:                          Options specific to Callisto1 or Callisto2 features (e.g., link up/down behavior, prefill).


"""
Base.@kwdef mutable struct Cal1Opts
    graph = Topology.Graph(["triangle"]; bidirectional = true)
    latency = 5000
    offset = :beta0
    tmax = 1e9
    kp = 2e-8
    ki = 1e-15
    poll_period = 100000
    theta0 = 0.1
    errors = :random
    base_freq = 1.0
    stopper = nothing
    check_freq_is_positive = true
    enable_logging = true
    controller = :pi
    # options below only used for callisto 1
    beta0 = 50
    gears = 1
    control_delay = 1000
    # options below written to by simulator
    model = nothing
    computed = nothing
end

# CalOpts should be such that at any time before running callisto,
# you can change any of the entries of CalOpts. That is, the constructor
# should not enforce constraints between the fields. That is why
# we do not instantiate the errors with a random vector of length
# n, since if subsequently the graph is changed then the struct
# would become inconsistent.
#
# Another viewpoint. Some parameters are set by the user
# and some are set by the simulator.
# No parameters should be set by both
#
#
# Also, dependencies of Callisto, like BittideTools, cannot access
# Cal1Opts.
#



function CalOpts(; topology = nothing, kw...)
    c = Cal1Opts(; allowed_kws(Cal1Opts, kw)...)
    if !isnothing(topology)
        c.graph = Topology.Graph(topology; bidirectional = true)
    end
    compute_parameters(c)
    return c
end


function make_ugn(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst)
    return beta0 - floor(gear*(theta0_at_src - latency * wm2_at_src)) + floor(gear*theta0_at_dst)
end

function compute_parameters(c::Cal1Opts)
    n = c.graph.n
    m = c.graph.m
    if c.errors == :random
        errors = make_frequencies(1, c.graph.n)
    else
        errors = c.errors
    end
    latency = [ati(c.latency, e) for e = 1:m]

    # beta0 is used to provide the measurement offsets
    # and to make the ugns
    if c.offset == :beta0
        offset = c.beta0
    else
        offset = c.offset
    end

    beta0 = c.beta0

    theta0 = [ati(c.theta0, i) for i = 1:n]
    gears = [ati(c.gears, e) for e=1:m]
    beta0 = [ati(beta0, e) for e=1:m]
    offset = [ati(offset, e) for e=1:m]

    if c.enable_logging
        # todo add recommended log length
        slog = log_init(m)
    else
        slog = nothing
    end
    if c.check_freq_is_positive
        checker = Checker()
    else
        checker = nothing
    end

    local_offsets = [inport_indexed_from_edge_indexed(c.graph, offset, i) for i=1:n]

    output_systems = [OutputSystem(local_offsets[i]) for i=1:n]

    if c.controller == :pi
        controller = [composed_controller(i, c.graph, c.kp,
            c.ki, c.poll_period, c.base_freq, output_systems[i]) for i=1:n]
    else
        controller = c.controller
    end
    epoch = -10*maximum(latency)
    incoming = c.graph.incoming_edges

    wm2 = [initial_error(errors[i]) for i=1:n]
    wm1 = [initial_error(errors[i]) for i=1:n]
    wmin = 0.9 * c.base_freq

    ugn =  [make_ugn(beta0[e], gears[e], latency[e], theta0[c.graph.edges[e].src],
                        wm2[c.graph.edges[e].src], theta0[c.graph.edges[e].dst]) for e=1:m]

    links = [Link(e, c.graph.edges[e].src, c.graph.edges[e].dst,
                 ugn[e], latency[e], gears[e], beta0[e], offset[e]) for e in 1:m]


    c.model = Model(;
        slog,
        tmax = c.tmax,
        epoch,
        num_nodes = n,
        links, incoming,
        p = c.poll_period,
        d = c.control_delay,
        errors,
        theta0,
        wm2,
        wm1,
        wmin,
        stopper = c.stopper,
        controller,
        checker)
    c.computed = ComputedParams(output_systems, offset, local_offsets, latency)
end



end
