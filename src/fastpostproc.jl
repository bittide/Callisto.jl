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


"""
    The module FastPostProc provides:

        `callistox`  which calls fast_post_process to generate the FastPost structure

"""
module FastPostProc


# The naming is all wrong
#
# This is the main post processing code for Callisto. So it should
# be in a module called Post.
#
# The function `callistox` is now used everywhere, so don't change this.
# But the FastPost struct should be something like PostProcessedSimData
# and the  fast_post_process function is really a constructor for ot.
#
# And we should remove all unused functions.
#



import ..Sim: callisto
import ..SimCore
import LogData: get_field, num_records
import Piecewise: Samples, differentiate, before, after

export  FastPost, fast_post_process, quick_post, callistox

Base.@kwdef mutable struct FastPost
    c
    simlog
    theta
    freq = nothing
    mocc = nothing
end

FastPost(c, simlog, theta) = FastPost(; c, simlog, theta)

#
# TODO: Move truncate to the Piecewise package.
#
function truncate(p, n)
    T = typeof(p)
    return T(p.x[1:end-n], p.y[1:end-n])
end

#
# TODO: determined whether truncate is needed here.
#
function fast_post_process(c, simlog, theta)
    fpost = FastPost(c, simlog, theta)
    # fpost.freq = truncate.(differentiate.(fpost.theta), 1)
    fpost.freq = after.(before.(differentiate.(fpost.theta), c.tmax), 0)
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
# Or just remove this function, as it is unused.
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
# TODO: remove this function, as it is unused.
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
