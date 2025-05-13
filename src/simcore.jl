

module SimCore

export Error, beta, gamma, initial_error, errorat, controllerfunctions

import ..Piecewise: PiecewiseLinear, definite_integral

# 
# abstract type Controller
# end
# 
# init(c::Controller) = return
# getlog(c::Controller) = 0.0
# step(c::Controller, i, theta_at_iofc, physical_time_at_iom, indegree, occ_measurements) = 0.0
# 
# mutable struct PIController <: Controller
#     kp::Float64
#     ki::Float64
#     xi::Float64
#     poll_period::Int64
#     base_freq::Float64
# end
# 
# function init(ctrl::PIController)
#     ctrl.xi = 0
# end
# 
# getlog(ctrl::PIController) = ctrl.xi
# 
# function step(ctrl::PIController, i::Int64,
#                            theta_at_instant_of_freq_change::Float64,
#                            physical_time_at_instant_of_measurement::Float64,
#                            indegree::Int64, occ_measurements::Vector{Float64})
#     #    r = sum(occ_measurements[1:indegree])
#     r = sum_some(occ_measurements, indegree)
#     next_state = ctrl.xi + ctrl.poll_period * r/ctrl.base_freq
#     correction = ctrl.kp*r + ctrl.ki * next_state
#     ctrl.xi = next_state
#     return correction
# end
# 







# uncontrolled frequency
abstract type Error end

struct ConstError <: Error
    const_error::Float64
end

struct PwlError <: Error
    pwl::PiecewiseLinear
end


initial_error(e::Error) = e(0)
initial_error(e::Number) = e
initial_error(e::Vector{T}) where {T <: Error} = e(0)
initial_error(e::Vector{Float64}) = e

errorat(e::Error, t) = e(t)
errorat(e::Number, t) = e
errorat(e::Vector{T}, t) where {T <: Number}  = e


Error(x::Number) = x
Error(p::PiecewiseLinear) = PwlError(p)

(e::ConstError)(t) = e.const_error
(e::PwlError)(t) = e.pwl(t)

(e::Vector{T})(t) where {T <: Error} = [a(t) for a in e]

Base.:+(a::ConstError, b::PiecewiseLinear) = Error(a.const_error + b)

# a function which takes
# node i, local time interval p, correction c, and wall-clock time s
# and returns a wall-clock time interval ds
#
# so that int_s^{s+ds} (c + \omega_i(t)) dt = p
#
# If omega_i is constant, then
#
#   ds*(c + omega_i) = p
#
# that is, ds = p / (c + omega_i)
#
#local_to_realtime(e::ConstError, p, c, s, wmin) = p / (c + e.const_error)
local_to_realtime(e::Number, p, c, s, wmin) = p / (c + e)
local_to_realtime(e::ConstError, p, c, s, wmin) = local_to_realtime(e.const_error, p, c, s, wmin)

# assume f(a) < 0 and f(b) > 0
function bisection(f, a, b, tol = 1e-5)
    while true
        c = (a + b)/2
        if abs(a-b) < tol
            return c
        end
        y = f(c)
        if y == 0
            return c
        end
        if y < 0
            a = c
        else
            b = c
        end
    end
end

function local_to_realtime(e::PwlError, p, c, s, wmin)
    #   intfreq(s, ds, c) = definite_integral(c + e.pwl, s, s + ds)
    #   dt = bisection(ds -> intfreq(s, ds, c) - p,  0, tmax)
    intfreq(s1, s2, c) = definite_integral(c + e.pwl, s1, s2)
    dt = bisection(s2 -> intfreq(s, s2, c) - p,  s, s + p/wmin) - s
    return dt
end


#beta((;ugn, src, dst, latency, gear), t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))

#function beta((;ugn, src, dst, latency, gear), t, theta)
#    println((;t, src, dst))
#    return ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))
#end

#gamma((;src, dst, latency, gear), t, theta) = floor(gear*theta[src](t)) - floor(gear*theta[src](t-latency))


##############################################################################




end
