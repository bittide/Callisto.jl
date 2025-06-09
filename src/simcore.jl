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

export Link

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

Base.@kwdef mutable struct C1Opts
    gears = 1
    beta0 = :offset
    control_delay = 1000
end



end
