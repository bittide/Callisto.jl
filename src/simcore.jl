

module SimCore

export Error, beta, gamma

import ..Piecewise: PiecewiseLinear


# uncontrolled frequency
abstract type Error end

struct ConstError <: Error
    const_error::Float64
end

struct PwlError <: Error
    pwl::PiecewiseLinear
end

Error(x::Number) = ConstError(x)
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
local_to_realtime(e::ConstError, p, c, s, wmin) = p / (c + e.const_error)
    
function local_to_realtime(e::PwlError, p, c, s, wmin)
    #   intfreq(s, ds, c) = definite_integral(c + e.pwl, s, s + ds)
    #   dt = bisection(ds -> intfreq(s, ds, c) - p,  0, tmax)
    intfreq(s1, s2, c) = definite_integral(c + e.pwl, s1, s2)
    dt = bisection(s2 -> intfreq(s, s2, c) - p,  s, s + p/wmin) - s
    return dt
end


beta((;ugn, src, dst, latency, gear), t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))

gamma((;src, dst, latency, gear), t, theta) = floor(gear*theta[src](t)) - floor(gear*theta[src](t-latency))






end
