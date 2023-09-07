
module Piecewise


export PiecewiseConstant, PiecewiseLinear, Series, Samples,
    integer_crossings, evaluate, differentiate,  delay, integrate, after,
    before, IntegersBetween, stairs, tuples

export crossing, invert, discontinuities, square, l2norm,
    integer_crossings_in_interval, definite_integral

#    steady_state, remove_steady_state, 
#    limits, xmax, 


abstract type Series end


"""
    PiecewiseConstant(x,y)

Struct for right continuous piecewise constant functions.

   f(t) = y[i] for t in [ x[i], x[i+1] )

These functions are defined for t in [x[1], x[end])
They do NOT include the endpoint. 
"""
struct PiecewiseConstant <: Series
    x
    y
    function PiecewiseConstant(x,y)
        # enforce invariants
        @assert is_increasing(x) "Error: in PiecewiseConstant, x must be increasing."
        @assert length(x) == length(y) + 1  "Error: in PiecewiseConstant, length(x) must equal length(y) + 1"
        return new(x,y)
    end
end





# piecewise linear function
#
#   theta(times[i]) = phase[i]   
#
# linear interpolation between points
#
# defined for all t in [times[1], times[end]]
# They DO include the endpoint.
# Note length(x) = length(y)
struct PiecewiseLinear <: Series
    x
    y
    function PiecewiseLinear(x,y)
        @assert is_increasing(x) "Error: in PiecewiseLinear, x must be increasing."
        @assert length(x) == length(y) "Error: in PiecewiseLinear, length(x) must equal length(y)"
        return new(x,y)
    end
end






# neither pwc nor pwl
struct Samples <: Series
    x
    y
    function Samples(x,y)
        @assert is_increasing(x) "Error: in Samples, x must be increasing."
        @assert length(x) == length(y) "Error: in Samples, length(x) must equal length(y)"
        return new(x,y)
    end
end


##############################################################################

# allow undef
#function uapply(p::Array, f)
#    out = Any[]
#    for i=1:length(p)
#        if isassigned(p,i)
#            push!(out, f(p[i]))
#        end
#    end
#    return out
#end


# Base.maximum(p::Series) = maximum(p.y)
# Base.minimum(p::Series) = minimum(p.y)

# 
# 
# function limits(p::Series) 
#     xmin = minimum(p.x)
#     xmax = maximum(p.x)
#     ymin = minimum(p.y)
#     ymax = maximum(p.y)
#     return xmin, xmax, ymin, ymax
# end

#function limits(p::Array{T}) where {T <: Series}
#    xmin = minimum(minimum.(uapply(p, a -> a.x)))
#    xmax = maximum(maximum.(uapply(p, a -> a.x)))
#    ymin = minimum(minimum.(uapply(p, a -> a.y)))
#    ymax = maximum(maximum.(uapply(p, a -> a.y)))
#    return xmin, xmax, ymin, ymax
#end

# 
# # p is an array 
# function limits(p::Array)
#     lims = Any[]
#     for ind in eachindex(skipmissing(p))
#         push!(lims, limits(p[ind]))
#     end
#     xmin = minimum(a[1] for a in lims)
#     xmax = maximum(a[2] for a in lims)
#     ymin = minimum(a[3] for a in lims)
#     ymax = maximum(a[4] for a in lims)
#     return xmin, xmax, ymin, ymax
# end
# 


"""
    tuples(p::PiecewiseLinear)

Return a list of tuples for plotting.
"""
tuples(p::PiecewiseLinear) = collect(zip(p.x, p.y))


"""
    tuples(p::PiecewiseLinear)

Return a list of tuples for plotting.
"""
tuples(p::Samples) = collect(zip(p.x, p.y))


tuples(p::Missing) = missing

# xmax(p::Series) = p.x[end]

# 
# function Base.push!(p::Series, a, b)
#     push!(p.x, a)
#    push!(p.y, b)
# end



##############################################################################



"""

    inv_interpolate(x1, y1, x2, y2, c)

Invert a linear interpolant f:x -> y between (x1,y1) and (x2,y2)
to find x such that f(x) = c

"""
inv_interpolate(x1,y1,x2,y2,c) = ((y2-c)*x1 - (y1-c)*x2)/(y2-y1)



"""

    interpolate(x1, y1, x2, y2, d)

Evaluate a linear interpolant f:x -> y between (x1,y1) and (x2,y2)
returning f(d).
"""
interpolate(x1,y1,x2,y2,d) = ((x2-d)*y1 - (x1-d)*y2)/(x2-x1)


"""

    interpolate(x::Vector, y::Vector, t)

Let f:x->y be the linear interpolant such that f(x[i]) = y[i].
Returns f(t).
"""
function interpolate(x::Vector, y::Vector, t)
    # this function returns the largest i such that x[i] <= t
    i = searchsortedlast(x, t)
    if x[i] == t
        return y[i]
    end
    @assert i != length(x)
    return interpolate(x[i], y[i], x[i+1], y[i+1], t)
end


evaluate(p::PiecewiseLinear, t) = interpolate(p.x, p.y, t)
invert(p::PiecewiseLinear, w) = interpolate(p.y, p.x, w)
(p::PiecewiseLinear)(t) = evaluate(p, t)


"""

    evaluate_with_lower_bound(p::PiecewiseLinear, t,  i0)

Given i0 such that p.x[i0] <= t, evaluate y, p(t).
Return y, i, where i is the largest integer such that p.x[i] <= t
"""
function  evaluate_with_lower_bound(p::PiecewiseLinear, t,  i0)
    i = findnext(a -> a >= t, p.x, i0)
    if p.x[i] == t
        return p.y[i], i
    end
    return interpolate(p.x[i-1], p.y[i-1], p.x[i], p.y[i], t), i 
end



"""
    tuples(p::PiecewiseConstant)

Return list of tuples giving points for plotting.
"""
function tuples(p::PiecewiseConstant)
    points = Vector{Tuple{Float64, Float64}}()
    for i=1:length(p.y)
        push!(points, (p.x[i],   p.y[i]))
        push!(points, (p.x[i+1], p.y[i]))
    end
    return points
end


#function xy(p::PiecewiseConstant)
#    x = Float64[]
#    y = Float64[]
#    function pt(a,b)
#        push!(x, a)
#        push!(y, b)
#    end
#    for i=1:length(p.y)
#        pt(p.x[i],   p.y[i])
#        pt(p.x[i+1], p.y[i])
#    end
#    return x, y
#end

# the points x,y where f(x)=y but f is discontinuous
# function discontinuities(p::PiecewiseConstant)
#    return p.x[1:end-1], p.y
# end

    
##############################################################################

# find integers i in the interval (xmin, xmax]
# This means in particular there are no integers in the interval
# when xmin == xmax, even if xmin and xmax are integral

"""
    IntegersBetween(a,b)

Iterator that returns the integers in the interval (a,b].
Note that (a,a] is always empty.
"""
struct IntegersBetween
    xmin
    xmax
end

Base.length(ib::IntegersBetween) = Int64(floor(ib.xmax) - floor(ib.xmin))
Base.zero(s::Series) = 0

function Base.iterate(ib::IntegersBetween)
    y = Int64(ceil(ib.xmin))
    if y == ib.xmin
        y = y + 1
    end
    if y <= ib.xmax
        return y, y
    end
    return nothing
end

function Base.iterate(ib::IntegersBetween, state)
    y = state + 1
    if y <= ib.xmax
        return y,y
    end
    return nothing
end


##############################################################################

"""
    integer_crossings(p::PiecewiseLinear)

Return integer values of p(t) where t in (xmin, xmax].
Assumes p is strictly increasing on that interval.
"""
function integer_crossings(p::PiecewiseLinear)
    m = length(IntegersBetween(p.y[1], p.y[end]))
    x = Array{Float64,1}(undef,m)
    y = Array{Float64,1}(undef,m)
    i = 1
    j = 1
    while i <= m
        for c in IntegersBetween(p.y[j], p.y[j+1])
            x[i] = inv_interpolate(p.x[j], p.y[j], p.x[j+1], p.y[j+1], c)
            y[i] = c
            i += 1
        end
        j += 1
    end
    return x, y
end



"""
    integer_crossings_in_interval(p::PiecewiseLinear, ymin, ymax)

Return (t,y) where y = p(t) and y in (ymin, ymax].
Assumes p is globally increasing.
"""
function integer_crossings_in_interval(p::PiecewiseLinear, ymin, ymax)
    # question only makes sense under these conditions
    @assert ymin >= p.y[1]
    @assert ymax <= p.y[end]

    m = length(IntegersBetween(ymin, ymax))
    x = Array{Float64,1}(undef,m)
    y = Array{Float64,1}(undef,m)
    i = 1
    # j = the index of the first value in p.y greater than or equal to ymin
    j = searchsortedfirst(p.y, ymin)
    # if j == 1 then p.y[1] == ymin so we don't need to look at
    # the interval (p.y[j-1], p.y[j]) because LHS of interval (ymin,ymax] is open
    #
    # otherwise we can have p.y[j] > ymin and so need to look at previous interval
    j = max(j - 1, 1)
    while i <= m
        lb = max(p.y[j], ymin)
        ub = min(p.y[j+1], ymax)
        for c in IntegersBetween(lb,ub)
            x[i] = inv_interpolate(p.x[j], p.y[j], p.x[j+1], p.y[j+1], c)
            y[i] = c
            i += 1
        end
        j += 1
    end
    # useless check that should always succeed unless there is a bug
    # so should be replaced by a test
    @assert i == length(x)+1
    return x, y
end


##############################################################################


"""
    floor(p::PiecewiseLinear)

Return a PiecewiseConstant function equal to the floor of p.
Assumes p is increasing. Note that floor is right continuous, so this can
return a PiecewiseConstant object.
"""
function Base.floor(p::PiecewiseLinear)
    # can have case where p = PiecewiseLinear([0, 100], [0.0, 0.0])
    # then xi = [], yi = []
    xi, yi = integer_crossings(p)
    pushfirst!(xi, p.x[1])
    pushfirst!(yi, yi[1]-1)
    if xi[end] < p.x[end]
        push!(xi, p.x[end])
    else
        yi = yi[1:end-1]
    end
    return PiecewiseConstant(xi, yi)
end

"""
    delay(p::PiecewiseLinear, tau)

Return q such that q(t) = p(t - tau)
"""
function delay(p::PiecewiseLinear, tau)
    x1 = p.x .+ tau
    y1 = copy(p.y)
    return PiecewiseLinear(x1, y1)
end


Base.:*(a::Number, pwc::PiecewiseConstant) = PiecewiseConstant(pwc.x, pwc.y*a)

function Base.:+(p::PiecewiseConstant, q::PiecewiseConstant)
    xmin = max(p.x[1], q.x[1])
    xmax = min(p.x[end], q.x[end])
    xn = unique(sort([p.x; q.x]))
    xn = [a for a in xn if xmin <= a <= xmax]
    n = length(xn)
    yn = zeros(n)
    pind = 1
    qind = 1
    for i=1:n-1
        pv, pind = evaluate_with_lower_bound(p, xn[i], pind)
        qv, qind = evaluate_with_lower_bound(q, xn[i], qind)
        yn[i] = pv + qv
    end
    return PiecewiseConstant(xn,yn[1:end-1])
end

function Base.:+(p::PiecewiseConstant, a::Number)
    xn = copy(p.x)
    yn = p.y .+ a
    return PiecewiseConstant(xn,yn)
end
Base.:+(a::Number, p::PiecewiseConstant) = p + a
Base.:-(p::PiecewiseConstant, q::PiecewiseConstant) =  p  + ((-1) * q)
Base.:-(p::PiecewiseConstant, a::Number) =  p  + ((-1) * a)
Base.:+(p::Samples, q::Samples)= Samples(p.x, p.y .+ q.y)
Base.:+(p::Samples, a::Number) = Samples(p.x, p.y .+ a)
Base.:+(a::Number, p::Samples) = p+a
Base.:-(p::Samples, a::Number) = Samples(p.x, p.y .- a)
Base.:-(a::Number, p::Samples) = Samples(p.x, a .- p.y)


##############################################################################

"""
    after(p::PiecewiseLinear, t)

Given p defined on interval [a,b].
Return q equal to p restricted to the interval [ max(a,t), b]
Requires t < b.
"""
function after(p::PiecewiseLinear, t)
    @assert t < p.x[end]  "Error: attempted to truncate p at too late a time"
    i = findfirst(a -> a >= t, p.x)
    if p.x[i] == t
        return PiecewiseLinear(p.x[i:end], p.y[i:end])
    end
    if i == 1
        return p
    end
    y = evaluate(p, t)
    return PiecewiseLinear([t ; p.x[i:end]], [y; p.y[i:end]])
end

"""
    before(p::PiecewiseLinear, t)

Given p defined on interval [a,b].
Return q equal to p restricted to the interval [ a, min(b,t) ]
Requires t > a.
"""
function before(p::PiecewiseLinear, t)
    @assert t > p.x[1] "Error: attempted to truncate p at too early a time"
    i = findlast(a -> a <= t, p.x)
    if p.x[i] == t
        return PiecewiseLinear(p.x[1:i], p.y[1:i])
    end
    if i == length(p.x)
        return p
    end
    y = evaluate(p, t)
    return PiecewiseLinear([p.x[1:i]; t], [p.y[1:i]; y])
end

"""
    integrate(p::PiecewiseConstant)

Given p defined on [a,b] returns q defined on [a,b]
defined by q(t) = int_{a}^{t} p(t) dt.
q is a PiecewiseLinear object.
"""
function integrate(p::PiecewiseConstant)
    x = p.x
    y = p.y
    n = length(x)
    z = Array{Float64,1}(undef, n)
    z[1] = 0
    for i=1:n-1
        z[i+1] = z[i] + y[i]*(x[i+1]-x[i])
    end
    xn = copy(p.x)
    return PiecewiseLinear(xn, z)
end

trapezoidal_integral(x1, y1, x2, y2) = (x2 - x1) * (y1 + y2) / 2

function integratex(x, y)
    n = length(x)
    z = Array{Float64,1}(undef, n)
    z[1] = 0
    for i=1:n-1
        #        dx = x[i+1] - x[i]
        #        z[i+1] = z[i] + (y[i] + y[i+1])*dx/2
        z[i+1] = z[i] + trapezoidal_integral(x[i], y[i], x[i+1], y[i+1])
    end
    xn = copy(x)
    return Samples(xn, z)
end


"""
    integrate(p::PiecewiseLinear)

Given p defined on [a,b] returns q defined on [a,b],
sampled at the same points as p. The function
q is quadratic, given by  q(t) = int_{a}^{t} p(t) dt.
"""
integrate(p::PiecewiseLinear) = integratex(p.x, p.y)


"""
    integrate(p::Samples)

Given p defined on [a,b] returns q defined on [a,b],
sampled at the same points as p. The function
q is given by  q(t) = nint_{a}^{t} p(t) dt,
computed using the trapezoidal approximation.
"""
integrate(p::Samples) = integratex(p.x, p.y)


"""
    differentiate(p::PiecewiseLinear)

Given p, return its derivative q, which is piecewise constant.
"""
function differentiate(p::PiecewiseLinear)
    x = p.x
    y = p.y
    n = length(x)
    z = Array{Float64,1}(undef, n-1)
    for i=1:n-1
        z[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
    end
    xn = copy(p.x)
    return PiecewiseConstant(xn, z)
end



"""

    evaluate_with_lower_bound(p::PiecewiseConstant, t,  i0)

Given i0 such that p.x[i0] <= t, evaluate y, p(t).
Return y, i, where i is the largest integer such that p.x[i] <= t
"""
function  evaluate_with_lower_bound(p::PiecewiseConstant, t,  i0)
    i = findnext(a -> a >= t, p.x, i0)
    if p.x[i] == t
        return p.y[i], i
    end
    return p.y[i-1], i
end

"""
    evaluate(p::PiecewiseConstant, t)

Return p(t). Note that by definition, p is defined
on an interval [a,b), where the right-hand end
of the interval is open. 
"""
function evaluate(p::PiecewiseConstant, t)
    i = findfirst(a -> a >= t, p.x)
    @assert !isnothing(i) "attempted to evaluate PiecewiseConstant(x) for x too large."
    @assert p.x[1] <= t   "attempted to evaluate PiecewiseConstant(x) for x too small."
    if i == length(p.x) && p.x[i] == t 
        error("attempted to evaluate PiecewiseConstant(x) for x at right hand endpoint.")
    end
    if p.x[i] == t
        return p.y[i]
    end
    if i == 1
        println("Error: attempted to evaluate p at too early a time")
        return 
    end
    return p.y[i-1]
end

# 
# function evaluate(p::PiecewiseConstant, t, allowrhend = false)
#     i = findfirst(a -> a >= t, p.x)
#     @assert !isnothing(i) "attempted to evaluate PiecewiseConstant(x) for x too large."
#     @assert p.x[1] <= t   "attempted to evaluate PiecewiseConstant(x) for x too small."
#     if i == length(p.x) && p.x[i] == t && !allowrhend
#         error("attempted to evaluate PiecewiseConstant(x) for x at right hand endpoint.")
#     end
#     if p.x[i] == t
#         return p.y[i]
#     end
#     if i == 1
#         println("Error: attempted to evaluate p at too early a time")
#         return 
#     end
#     return p.y[i-1]
# end
# 
# 

(p::PiecewiseConstant)(t) =  evaluate(p, t)



"""
    delay(p::PiecewiseConstant, tau)

Return q such that q(t) = p(t - tau)
"""
function delay(p::PiecewiseConstant, tau)
    xn = p.x .+ tau
    yn = copy(p.y)
    return PiecewiseConstant(xn, yn)
end


##############################################################################


squarex(p::T) where {T<:Series} =  T(p.x, p.y .* p.y)
square(p::Samples) =  squarex(p)
square(p::PiecewiseConstant) =  squarex(p)

l2norm(p::Series) = integrate(square(p)).y[end]
l2norm(p::Array{T,1}) where {T<:Series} = sum(l2norm.(p))


# steady_state(p::Series) =  p.y[end]
# steady_state(p::Array{T,1}) where {T<:Series} = [steady_state(p[i]) for i=1:length(p)]
# remove_steady_state(p::Series) = p - steady_state(p)
# relative_l2norm(p::Series) = l2norm(remove_steady_state(p))
# export l2norm, steady_state, remove_steady_state, relative_l2norm
# l2norm(p::Series, offset) = l2norm(p .- offset)
    
Base.:*(a::Number, pwc::PiecewiseLinear) = PiecewiseLinear(pwc.x, pwc.y*a)
Base.:*(a::Number, p::Samples) = Samples(p.x, p.y*a)
Base.:+(p::PiecewiseLinear, a::Number) = PiecewiseLinear(copy(p.x), p.y .+ a)
Base.:+(a::Number, p::PiecewiseLinear) = p + a
Base.:-(p::PiecewiseLinear, q::PiecewiseLinear) =  p  + ((-1) * q)
Base.:-(p::PiecewiseLinear, a::Number) =  p  + ((-1) * a)


function Base.:+(p::PiecewiseLinear, q::PiecewiseLinear)
    xmin = max(p.x[1], q.x[1])
    xmax = min(p.x[end], q.x[end])
    xn = unique(sort([p.x; q.x]))
    xn = [a for a in xn if xmin <= a <= xmax]
    n = length(xn)
    yn = zeros(n)
    pind = 1
    qind = 1
    for i=1:n
        pv, pind = evaluate_with_lower_bound(p, xn[i], pind)
        qv, qind = evaluate_with_lower_bound(q, xn[i], qind)
        yn[i] = pv + qv
    end
    return PiecewiseLinear(xn,yn)
end




##############################################################################


"""
    after(p::PiecewiseConstant, t)

Given p defined on interval [a,b].
Return q equal to p restricted to the interval [ max(a,t), b]
Requires t < b.
"""
function after(p::PiecewiseConstant, t)
    @assert t < p.x[end]  "Error: attempted to truncate p at too late a time"
    i = findfirst(a -> a >= t, p.x)
    if p.x[i] == t
        return PiecewiseConstant(p.x[i:end], p.y[i:end])
    end
    if i == 1
        return p
    end
    y = evaluate(p, t)
    return PiecewiseConstant([t ; p.x[i:end]], [y; p.y[i:end]])
end



"""
    before(p::PiecewiseConstant, t)

Given p defined on interval [a,b].
Return q equal to p restricted to the interval [ a, min(b,t) ]
Requires t > a.
"""
function before(p::PiecewiseConstant, t)
    @assert t > p.x[1] "Error: attempted to truncate p at too early a time"
    i = findlast(a -> a <= t, p.x)
    if p.x[i] == t
        return PiecewiseConstant(p.x[1:i], p.y[1:i-1])
    end
    if i == length(p.x)
        return p
    end
    y = evaluate(p, t)
    return PiecewiseConstant([p.x[1:i]; t], [p.y[1:i-1]; y])
end


"""
    before(p::Samples, t)

Given p defined on interval [a,b].
Return q equal to p restricted to the interval [ a, min(b,t) ]
Requires t > a.
"""
function before(p::Samples, t)
    @assert t > p.x[1] "Error: attempted to truncate p at too early a time"
    i = findlast(a -> a <= t, p.x)
    return Samples(p.x[1:i], p.y[1:i])
end



"""
    is_nondecreasing(x)

Return true if x[i-1] <= x[i] for all i=2...length(x).
"""
function is_nondecreasing(x)
    for i=2:length(x)
        if x[i-1] > x[i]
            return false
        end
    end
    return true
end
    

"""
    is_increasing(x)

Return true if x[i-1] < x[i] for all i=2...length(x).
"""
function is_increasing(x)
    for i=2:length(x)
        if x[i-1] >= x[i]
            return false
        end
    end
    return true
end


##############################################################################
# definite integral


"""
    definite_integral(p::PiecewiseLinear, x1, x2)

Return the integral of p over the interval [x1, x2].
"""
function definite_integral(p::PiecewiseLinear, x1, x2)
    @assert x1 >= p.x[1] "Left integration x1 = $(x1) limit is too small"
    @assert x2 < p.x[end] "Left integration limit is too large"
    if x2 == x1
        return 0
    end
    i0 = findlast(a -> a <= x1, p.x)
    ev(x) = evaluate_with_lower_bound(p, x, i0)[1]
    i = i0
#    if isnothing(i)
#        println("Left integration x1 = $(x1) limit is too small")
#        return nothing
#    end
#    if i == length(p.x)
#        println("Left integration limit is too large")
#        return nothing
#    end
    integral_so_far = 0.0
    x = x1
    while p.x[i+1] <= x2
        y = ev(x)
        integral_so_far += trapezoidal_integral(x, y, p.x[i+1], p.y[i+1])
        x = p.x[i+1]
        i += 1
    end
    if p.x[i] < x2
        y = ev(x)
        y2 = ev(x2)
        integral_so_far += trapezoidal_integral(x, y, x2, y2)
    end
    return integral_so_far
end


end
