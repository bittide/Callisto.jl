
module FastPostProc

import ..SimCore
import ..LogData:    num_records, get_field
import ..Piecewise:  Samples, differentiate
import ..Sim:  callisto

export  FastPost, fast_post_process, quick_post, callistox

Base.@kwdef mutable struct FastPost
    c
    simlog
    theta
    freq = nothing
    mocc = nothing
end

FastPost(c, simlog, theta) = FastPost(; c, simlog, theta)

function truncate(p, n)
    T = typeof(p)
    return T(p.x[1:end-n], p.y[1:end-n])
end


function fast_post_process(c, simlog, theta)
    fpost = FastPost(c, simlog, theta)
    fpost.freq = truncate.(differentiate.(fpost.theta), 1)
    fpost.mocc =  [get_mocc(fpost, e) for e=1:c.graph.m]
    return fpost
end


function get_mocc(fpost, e)
    tk = get_field(fpost.simlog[e], :tk)
    mocc = get_field(fpost.simlog[e], :mocc)
    return Samples(tk, mocc)
end


callistox(c) = fast_post_process(c, callisto(c)...)


##############################################################################
# For those who are extra hasty, use this.
# This is not approximate, but it only provides equidistant sampled values.


# TODO: move this to a generic function in Piecewise
# that returns derivative values at specific points.
#
#
# compute the sampled frequency by differentiating theta at tvals
#
function get_freq_samples(theta, tvals)
    freqs = []
    for t in tvals
        i = searchsortedlast(theta.x, t)
        if i == length(theta.x)
            i = i - 1
        end
        push!(freqs, (theta.y[i+1] - theta.y[i])/(theta.x[i+1] - theta.x[i]))
    end
    return freqs
end


#
# Sample the occupancy at N equidistant points between 0 and the
# max simulation time. Also sample the frequency, by differentiating
# theta, instead of using the logged value.
#
function quick_post(c, theta; N = 500)
    tmax = minimum(a.x[end] for a in theta)
    n = length(theta)
    tvals = collect(range(0, tmax, length = N))
    freq = Samples[]
    for i = 1:n
        fr = get_freq_samples(theta[i], tvals)
        push!(freq, Samples(tvals, fr))
    end
    occ = get_occ_samples(c, theta, tvals)
    qc = (freq = freq, occ = occ)
    return qc
end


end
