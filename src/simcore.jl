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
using StateSystems
using DataStructures

export Cal1Opts, CalOpts, compute_parameters, get_beta0,
    get_gears, get_latency, get_offsets,
    get_output_systems, get_ugns, Link, log_init, log_occupancy, Event, Model, Status

export Frame, FrameMover, make_framemover

const log_fields = (; tk=1, mocc=2)
log_init(m) =  [Log(log_fields) for i=1:m]
@inline function log_occupancy(slog::Vector{Log}, eid, t, mocc)
    LogData.log(slog[eid], log_fields.tk, t)
    LogData.log(slog[eid], log_fields.mocc, mocc)
    start_new_record(slog[eid])
end
@inline log_occupancy(slog::Nothing, eid, t, mocc) = return




struct Frame
    send_tick::UInt64
    data::UInt64
end

struct FrameMover
    fifos::Vector{DataStructures.CircularBuffer{Frame}}
end

function make_framemover(m)
    fifos = [DataStructures.CircularBuffer{Frame}(100) for e = 1:m]
    return FrameMover(fifos)
end




"""
    Model

Holds the computed parameters for a Callisto simulation.

Unlike `CalOpts`, which holds user-facing configuration, `Model` contains the
parameters directly used by the simulator. It is typically created by the
`compute_parameters` function from a `CalOpts` object.

# Fields
- `slog`: An optional vector of `Log` objects for recording buffer occupancy data for each edge.
- `tmax`: The maximum simulation time.
- `epoch`: A time before `t=0` used to initialize the simulation state, particularly the history of `theta`.
- `p`: The polling period in local clock ticks between controller updates.
- `d`: The control delay in local clock ticks between measurement and actuation.
- `errors`: A vector of clock frequency error models for each node.
- `theta0`: A vector of initial clock phases for each node at `t=0`.
- `wm2`: A vector of initial frequencies for epoch < t < 0.
- `wm1`: A vector of initial frequencies for 0 < t < s0, where s0 is the time when local clock = theta0 + d
- `wmin`: The minimum allowed frequency, used to prevent division by zero in calculations.
- `stopper`: An optional custom function or object that can terminate the simulation early.
- `controller`: A vector of controller objects, one for each node, that determine frequency corrections.
- `checker`: An optional `Checker` object to validate that clock frequencies remain positive.
- `ugns`: A vector of integer offsets for each link used in the buffer occupancy calculation.
- `graph`: The network topology as a `Topology.Graph` object.
- `latency`: A vector of latencies for each link in the graph.
- `gears`: A vector of gear ratios for each link, used in buffer occupancy calculations.
"""
Base.@kwdef mutable struct Model
    slog
    tmax
    epoch
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
    ugns
    graph
    latency
    gears
    framemover
    framedataprocessor
end


struct Event
    node::Int64
    time::Float64
end

Base.@kwdef mutable struct ComputedParams
    output_systems
    offset
    local_offsets
    latency
    ugns::Vector{Int64}
    gears::Vector{Float64}
    beta0::Vector{Int64}
end

"""
    Status

Holds the dynamic state of a Callisto simulation.

State variables of the simulation.
Initialized by `make_status` and then passed through the simulation loop,
where it is updated at each simulation event.

# Fields
- `model::Model`: The static parameters of the simulation model.
- `measurement::Vector{Measurement}`: A vector where `measurement[i]` holds the most
  recent measurement data for node `i`. This includes buffer occupancies for all
  incoming edges, the local clock value at the instance of control (`theta_at_ioc`),
  and the physical time of measurement (`physical_time_at_iom`).
- `theta`: A vector of `PiecewiseLinear` functions, where `theta[i]` represents the
  local clock `theta` of node `i` as a function of physical time.
- `guess1::Int64`: An integer value used as a heuristic to speed up the evaluation
  of `theta` when computing buffer occupancies. It provides an initial guess for
  look-behind time.
- `guess2::Int64`: Another integer guess, similar to `guess1`, used to accelerate
  the inversion of `theta` functions.
"""
Base.@kwdef mutable struct Status
    model::Model
    measurement::Vector{Measurement}
    theta
    guess1::Int64
    guess2::Int64
end

"""
    Cal1Opts structure
Holds configuration options for a Callisto simulation.
Fields
- `graph`:                      The network topology. Defaults to a 3-node bidirectional triangle.
- `latency`:                    Link latencies. Can be a scalar (applied to all links) or a vector. Defaults to 5000.
- `offset`:                     Initial buffer offsets/target occupancies. Can be a scalar or a vector, or `:beta0` to use `beta0`. Defaults to `:beta0`.
- `tmax`:                       Maximum simulation time. Defaults to 1e9.
- `kp`:                         Proportional gain for the PI controller. Defaults to 2e-8.
- `ki`:                         Integral gain for the PI controller. Defaults to 1e-15.
- `poll_period`:                Time between controller updates for each node. Defaults to 100000.
- `theta0`:                     Initial phase for each node. Can be a scalar or a vector. Defaults to 0.1.
- `errors`:                     Clock frequency errors for each node. Can be `:random` (generates random errors), a vector of numerical error values, or a vector of `SimCore.Error` objects. Defaults to `:random`.
- `base_freq`:                  Base frequency used by the controller for scaling the integral term. Defaults to 1.0.
- `stopper`:                    A custom stopping condition for the simulation. Defaults to `nothing` (no custom stopper).
- `check_freq_is_positive`:     If `true`, checks if node frequencies become negative during simulation. Defaults to `true`.
- `enable_logging`:             If `true`, enables logging of simulation data (e.g., buffer occupancies). Defaults to `true`.
- `controller`:                 Specifies the controller type for each node. Can be `:pi` (default PI controller) or a vector of custom controller objects. Defaults to `:pi`.
- `beta0`:                      The nominal value of beta, used by the UGN calculation. Can be a scalar or a vector. Defaults to 50.
- `gears`:                      The gear ratio between the nominal phase increment of connected nodes. Can be a scalar or a vector. Defaults to 1.
- `control_delay`:              The control delay in local clock ticks between measurement and actuation. Defaults to 1000.
- `model`:                      Field written to by `compute_parameters` and used by the simulator.
- `computed`:                   Field written to by `compute_parameters` but not used by the simulator.

CalOpts is such that at any time before running callisto,  you can change any of the entries of CalOpts.
That is, the constructor should not enforce constraints between the fields. That is why
we do not instantiate the errors with a random vector of length n, since if subsequently
the graph is changed then the struct would become inconsistent.

Most fields are parameters are to be set by the user. The `computed` and `model` fields are
set by the simulator or initializer.

Dependencies of Callisto, like BittideTools, cannot access Cal1Opts.
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
    framemover = nothing
    framedataprocessor = nothing
    # options below written to by simulator
    model = nothing
    computed = nothing
end

######################################################################################
# Functions for the user
get_ugns(c::Cal1Opts) = c.computed.ugns
get_latency(c::Cal1Opts) = c.computed.latency
get_gears(c::Cal1Opts) = c.computed.gears
get_output_systems(c::Cal1Opts) = c.computed.output_systems
get_offsets(c::Cal1Opts) = c.computed.offset
get_beta0(c::Cal1Opts) = c.computed.beta0
###########################################################################################

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

    # beta0 is used to provide the measurement offsets and to make the ugns
    if c.offset == :beta0
        offset = c.beta0
    else
        offset = c.offset
    end

    theta0 = [ati(c.theta0, i) for i = 1:n]
    gears = [ati(c.gears, e) for e=1:m]
    beta0 = [ati(c.beta0, e) for e=1:m]
    offset = [ati(offset, e) for e=1:m]

    if c.enable_logging
        # TODO: add recommended log length
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
        controller = [compose(PIStateSystem(c.kp, c.ki * c.poll_period / c.base_freq), output_systems[i]) for i=1:n]
    else
        controller = c.controller
    end
    epoch = -10*maximum(latency)

    wm2 = [initial_error(errors[i]) for i=1:n]
    wm1 = [initial_error(errors[i]) for i=1:n]
    wmin = 0.9 * c.base_freq

    ugns = [make_ugn(beta0[e], gears[e], latency[e], theta0[c.graph.edges[e].src],
                        wm2[c.graph.edges[e].src], theta0[c.graph.edges[e].dst]) for e=1:m]


    if c.framemover == :fifo
        framemover = make_framemover(c.graph.m)
    elseif isnothing(c.framemover)
        framemover = nothing
    else
        framemover = c.framemover
    end

    c.model = Model(;
        slog,
        tmax = c.tmax,
        epoch,
        p = c.poll_period,
        d = c.control_delay,
        errors,
        theta0,
        wm2,
        wm1,
        wmin,
        stopper = c.stopper,
        controller,
        checker,
        ugns,
        graph = c.graph,
        latency,
        gears,
        framemover,
        framedataprocessor = c.framedataprocessor)

    c.computed = ComputedParams(output_systems, offset, local_offsets, latency, ugns, gears, beta0)
end



end
