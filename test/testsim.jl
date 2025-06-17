module TestSim

using Piecewise
using Topology
using Callisto.Sim
using Callisto.SimCore
using BittideTools
using StateSystems

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

function make_ugn(beta0, gear, latency, theta0_at_src,  wm2_at_src, theta0_at_dst)
    return beta0 - floor(gear*(theta0_at_src - latency * wm2_at_src)) + floor(gear*theta0_at_dst)
end

function testcallisto(tmax,
                      graph,
                      poll_period,
                      kp,
                      ki,
                      errors,
                      base_freq)
    m = graph.m
    n = graph.n
    theta0 = fill(0.1, n)
    wm2 = copy(errors)
    wm1 = copy(errors)
    beta0 = fill(50, m)
    latency = 5000
    gear = 1
    offset = fill(50, m)

    function make_link(e)
        src = graph.edges[e].src
        dst = graph.edges[e].dst
        ugn = make_ugn(beta0[e], gear, latency, theta0[src], wm2[src], theta0[dst])
        return Link(e, src, dst, ugn, latency, gear, beta0[e], offset[e])
    end

    links = [make_link(e) for e = 1:m]
    wmin = 0.1
    epoch = -10*latency
    control_delay = 1000
    slog = nothing
    checker = nothing
    stopper = nothing

    local_offsets = [inport_indexed_from_edge_indexed(graph, offset, i) for i=1:n]
    output_systems = [OutputSystem(local_offsets[i]) for i=1:n]
    controller = [compose(PIStateSystem(c.kp, c.ki * c.poll_period / c.base_freq), output_systems[i]) for i=1:n]


    simlog, theta = callisto(slog, tmax, epoch, graph.n,
        links, graph.incoming_edges, poll_period, control_delay,
        errors, theta0, wm2, wm1, wmin, stopper,
        controller, checker)

    return theta
end


function main1()
    kp = 2e-8
    ki = 1e-15
    tmax = 200e6
    alpha = 1e-4
    f0 = 1
    freqs = [f0 + alpha, f0-alpha, f0]
    base_freq = f0
    graph = Topology.Graph(["triangle"]; bidirectional = true)
    poll_period = 100000
    theta =  testcallisto(tmax, graph, poll_period, kp, ki, freqs, base_freq)
    end_time_in_range = all([a.x[end] > tmax for a in theta])
    ns = sum(offsetnormsquared(q, f0) for q in theta)
    # resistance distance test
    res = 2/3
    expected_ns = alpha^2*res/(2*kp)
    norm_of_freq_in_range = 0.9 <= ns/expected_ns <= 1.1
    return norm_of_freq_in_range && end_time_in_range
end

function main2()
   return 1 == 1
end

function main3()
    return 1 == 1
end



end
