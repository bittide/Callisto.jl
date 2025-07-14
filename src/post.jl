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
   The module Post provides the function `focused_callisto_info` which
   provides the exact buffer occupancy as a function of time.
"""
module Post

using Piecewise

export focused_callisto_info, get_occupancy, get_occupancies, get_sampled_occupancies

import Piecewise: after, before, delay
import ..SimCore: get_latency, get_ugns, get_gears, Cal1Opts

#
# these functions depend on the sim results (i.e, theta)
# and so they are here instead of in simcore.
#
function get_occupancies(c::Cal1Opts, theta, t::Number)
    latency = get_latency(c)
    gears = get_gears(c)
    ugns = get_ugns(c)
    return [get_occupancy(ugns[e], c.graph.edges[e], latency[e], gears[e], theta, t) for e = 1:c.graph.m]
end

function get_sampled_occupancies(c::Cal1Opts, theta, tlist::Vector)
    occupancies_plus_times = Array{Samples, 1}(undef, c.graph.m)
    for e = 1:c.graph.m
        occupancies = [get_occupancy(c, theta, t, e) for t in tlist]
        occupancies_plus_times[e] = Samples(tlist, occupancies)
    end
    return occupancies_plus_times
end

get_occupancy(ugn, edge, latency, gear, theta, t) = ugn + floor(gear*theta[edge.src](t - latency)) - floor(gear*theta[edge.dst](t))

function get_occupancy(c, theta, t, e)
    latency = get_latency(c)
    gears = get_gears(c)
    ugns = get_ugns(c)
    return get_occupancy(ugns[e], c.graph.edges[e], latency[e], gears[e], theta, t)
end



##############################################################################
# For those who want all the details, use this

# buffer at node i associated with link to node j
# find occupancy as a piecewise constant function
function beta_piecewise(edge, latency, ugn, gear, theta, tmin, tmax)
    thj1 = delay(theta[edge.src], latency)
    thj = before(after(thj1, tmin), tmax)
    thi = before(after(theta[edge.dst], tmin), tmax)
    p = floor(gear*thj) - floor(gear*thi)
    return p + ugn
end


# things that take a long time we only apply to a time range
function get_occ(c, theta, tmin, tmax)
    latency = get_latency(c)
    ugns = get_ugns(c)
    gears = get_gears(c)
    occ = [beta_piecewise(c.graph.edges[e], latency[e], ugns[e], gears[e], theta, tmin, tmax) for e = 1:c.graph.m]
    return occ
end


function focused_callisto_info(c, xc, tmin, tmax)
    occ = get_occ(c, xc.theta, tmin, tmax)
    fc = (occ = occ, )
    return fc
end



end

