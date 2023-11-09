

module SimCore

export Error, beta, gamma

import ..Piecewise: PiecewiseLinear, definite_integral



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


beta((;ugn, src, dst, latency, gear), t, theta) = ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))

#function beta((;ugn, src, dst, latency, gear), t, theta)
#    println((;t, src, dst))
#    return ugn + floor(gear*theta[src](t - latency)) - floor(gear*theta[dst](t))
#end

gamma((;src, dst, latency, gear), t, theta) = floor(gear*theta[src](t)) - floor(gear*theta[src](t-latency))






end
