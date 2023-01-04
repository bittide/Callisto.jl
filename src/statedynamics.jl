
# code for state-space systems
module StateDynamics



##############################################################################
# asbtract

abstract type StateSystem
end

##############################################################################
# PI controller

mutable struct PIStateSystem <: StateSystem
    x
    kp
    ki
end

function PIStateSystem(kp, ki)
    return PIStateSystem(0, kp, ki)
end

function next(K::PIStateSystem, u)
    next_state =  K.x + u 
    y =  K.kp * u +  K.ki * next_state
    K.x = next_state
    return y
end


##############################################################################
# P controller


mutable struct StaticStateSystem <: StateSystem
    f
end


next(K::StaticStateSystem, u) = K.f(u)

##############################################################################

mutable struct StaticPlusConstantStateSystem <: StateSystem
    f
    constant
end


next(K::StaticPlusConstantStateSystem, u) = K.f(u) + K.constant

function setconstant(K::StaticPlusConstantStateSystem, c)
    K.constant = c
end


##############################################################################
# Constant controller

mutable struct ConstantStateSystem <: StateSystem
    constant
end


next(K::ConstantStateSystem, u) =  K.constant

function setconstant(K::ConstantStateSystem, c)
    K.constant = c
end


##############################################################################
# switching wrapper

mutable struct SequentialStateSystem
    count        # time
    lastoutput
    subsystems   # list of controller objects
    endtimes     # times on which to switch
    switches     # list of functions called at switch events

end

function SequentialStateSystem(subsystems, endtimes, switches)
    return SequentialStateSystem(0, 0, subsystems, endtimes, switches)
end

# endtimes = [5, 10]
#   which_subsystem = 1  for t <= 5
#   which_subsystem = 2  for 6 <= t <= 10
which_subsystem(endtimes, t) = searchsortedfirst(endtimes, t)


# controller t takes state t and measurement, determines state t+1
function next(K::SequentialStateSystem, measurement)
    t = K.count 

    previous_subsystem_index = which_subsystem(K.endtimes, t-1)
    this_subsystem_index = which_subsystem(K.endtimes, t)
    if this_subsystem_index > previous_subsystem_index
        K.switches[previous_subsystem_index](K.lastoutput)
    end
    K.count += 1
    y = next(K.subsystems[this_subsystem_index], measurement)
    K.lastoutput = y
    return y
 end

##############################################################################
# system that switches off

function StoppingStateSystem(toff, K1)
    K2 = ConstantStateSystem(0)
    function sf(lastoutput)
        K2.constant = lastoutput
    end
    K = SequentialStateSystem([K1, K2], [toff], [sf])
end


########################################################################
# composition of controllers

mutable struct CompositionStateSystem <: StateSystem
    K1
    K2
end

function next(K::CompositionStateSystem, u)
    y2 = next(K.K2, u)
    y = next(K.K1, y2)
    return y
end


Base.:*(K1::StateSystem, K2::StateSystem) = CompositionStateSystem(K1, K2)


########################################################################
# specific to bittide
# takes a measurement input
#

struct AverageOccupancyStateSystem <: StateSystem end

function next(K::AverageOccupancyStateSystem, measurement)
    if length(measurement) == 0
        # can have no measurements if
        # there are no incoming edges at a node
        # and hence no elastic buffers
        return 0
    end
    y = sum(a[2] for a in measurement)
    return y
end

########################################################################
# specific to bittide
# returns the occupancy of one elastic buffer


struct OneEdgeOccupancyStateSystem <: StateSystem
    edgeid
end

function next(K::OneEdgeOccupancyStateSystem, measurement)
    i = findfirst(x -> x[1] == K.edgeid,  measurement)
    if isnothing(i)
        println("Error: cannot find edge $(K.edgeid)")
        return 0
    end
    return measurement[i][2]
end


##############################################################################
# automatic reset controller, a variant of PI


mutable struct AutomaticResetStateSystem <: StateSystem
    x
    b
    alpha
end

function AutomaticResetStateSystem(b, alpha)
    return AutomaticResetStateSystem(0, b, alpha)
end

function next(K::AutomaticResetStateSystem, u)
    y =  K.x + K.b * u 
    K.x  = (1 - K.alpha) * K.x   +   K.alpha * y 
    return y
end

automaticresetcontroller(c) = AutomaticResetStateSystem(c.kp, c.ki * c.poll_period / c.kp) * AverageOccupancyStateSystem()



##############################################################################
# convenience functions

picontroller(c) = PIStateSystem(c.kp, c.ki * c.poll_period / c.base_freq) * AverageOccupancyStateSystem()
pcontroller(c) = StaticStateSystem(u -> c.kp*u) * AverageOccupancyStateSystem()


##############################################################################
# controller functions

getopts(; kw...) = kwtuple(kw)
kwtuple(opts) = NamedTuple{keys(opts)}(values(opts))

function controllerfunctions(controllers)
    controller_init = () -> controllers
    controller_next(i, Klist, measurements) =  (Klist, next(Klist[i], measurements))
    controller_log(i, Klist) = 0
    return getopts(; controller_init, controller_next, controller_log)
end




=======

end
