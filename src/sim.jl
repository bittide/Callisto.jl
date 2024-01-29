


module Sim

using DataStructures

export callisto

import ..Piecewise:  PiecewiseLinear, invert
import ..LogData:    Log, llog, llognext
import ..SimCore:    beta, local_to_realtime

const USE_HEAP = true

# also needed in opts.jl
# This is useful for very large simulations
# where we give up flexibility and logging for speed
const USE_SUM_MEASUREMENT = false
##############################################################################

function clog(slog, sk, nid, freq)
    llog(slog, 1, 1)
    llog(slog, 3, sk)
    llog(slog, 4, nid)
    llog(slog, 6, freq)
    llognext(slog)
end

function clog(slog, k, sk, nid, freq)
    llog(slog, 1, 1)
    llog(slog, 2, k)
    llog(slog, 3, sk)
    llog(slog, 4, nid)
    llog(slog, 6, freq)
    llognext(slog)
end

##############################################################################




meas(link, t, theta, betafn) = (betafn(link, t, theta) - link.offset)/link.gear

##############################################################################

function Gx(i, theta, s, theta0, p, links, incoming, d, betafn, slog)
    #
    # This means in particular that sampling always starts
    # off roughly in phase across all nodes. Perhaps we should have
    # an option where we can test different initial phasing.
    #
    k = round((theta[i].y[end] - d - theta0[i])/p)

    # following two approaches are both fine
    # t = invert(theta[i], theta[i].y[end] - d)
    t = invert(theta[i], theta0[i] + k*p)
    llog(slog, 2, k)
    llog(slog, 5, t)
    if USE_SUM_MEASUREMENT
        y = 0.0
        for e in incoming[i]
            y += meas(links[e], t, theta, betafn)
        end
    else
        y = [ (e, meas(links[e], t, theta, betafn), t, theta) for e in incoming[i]]
    end
    return y
end

function sum_or_missing(x)
    if length(x) == 0
        return -100
    end
    return sum(x)
end

# note controller cannot exactly implement PI, since it does not know t
function Cx(i, y, controller_states, controller_next, controller_log, slog)
    if USE_SUM_MEASUREMENT
        llog(slog, 7, y)
    else
        llog(slog, 7, sum_or_missing(a[2] for a in y))
    end
        
    llog(slog, 8, controller_log(i, controller_states[i]))
    controller_states[i], corr = controller_next(i, controller_states[i], y)
    return corr
end

function extend(p::PiecewiseLinear, dx, dy)
    push!(p.x, p.x[end] + dx)
    push!(p.y, p.y[end] + dy)
end

function Fx(i, theta, ds, p)
    # p/(c+e) is the time interval
    # between s_k and s_{k+1}
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


function callisto(tmax, epoch, num_nodes, links, incoming, p, d, errors,
                  theta0, freq_m2, freq_m1, wmin, betafn,
                  controller_init, controller_next, controller_log,
                  special, stopper, sim_is_done)

        
    @assert d < p
    if tmax/p < 10
        println("Warning; tmax < 10 * poll_period")
    end
    
    header = ["ltype", "k", "sk", "nid", "tk", "freq", "meas", "xi"]
    slog = Log(8,header)

    initial_state(i) = initial_statex(i, epoch, freq_m2[i], freq_m1[i], theta0[i], d, slog)
    G(i, theta, s)   = Gx(i, theta, s, theta0, p, links, incoming, d, betafn, slog)
    C(i, y)          = Cx(i, y, controller_states, controller_next, controller_log, slog)
    F(i, theta, ds)   = Fx(i, theta, ds, p)

    theta = [initial_state(i) for i=1:num_nodes]
    controller_states = [ controller_init() for i=1:num_nodes]
    s = 0.0
    if USE_HEAP
        h = MutableBinaryMinHeap{Float64}()
        for j=1:num_nodes
            i = push!(h, theta[j].x[end])
            @assert i == j
        end
    end
    done = false
    while s < tmax && !done
        if USE_HEAP
            s, i = top_with_handle(h)
        else
            s, i  = findmin([theta[j].x[end] for j=1:num_nodes])
        end
        y = G(i, theta, s)
        c = C(i, y)
        ds = local_to_realtime(errors[i], p, c, s, wmin)
        w = c + errors[i](s)
        clog(slog, theta[i].x[end], i, w)
        if c + errors[i](s) <= 0
            println("Frequency at node $i has become negative: omega = ", c + errors[i](s))
            println("time s = ", s)
            return slog, theta
        end
        #@assert c + errors[i](s) > 0
        F(i, theta, ds)
        if USE_HEAP
            update!(h, i, theta[i].x[end])
        end
        done = sim_is_done(stopper, s, i, w)
    end
    return slog, theta
end


function callisto(c)
    print("Running callisto: ")

    # callisto cares about only the incoming edges to each node
    @time simlog, theta = callisto(c.tmax,
                                   c.epoch,
                                   c.graph.n,
                                   c.links,
                                   c.graph.incoming_edges,
                                   c.poll_period,
                                   c.control_delay,
                                   c.errors,
                                   c.theta0,
                                   c.wm2,
                                   c.wm1,
                                   c.wmin,
                                   c.betafn,
                                   c.controller_init,
                                   c.controller_next,
                                   c.controller_log,
                                   c.special,
                                   c.stopper,
                                   c.sim_is_done)

    return (simlog = simlog, theta = theta,)
end










end
