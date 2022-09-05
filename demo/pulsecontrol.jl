
module Run

using Callisto
using PlotKit
const pk = PlotKit
using Random

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)
callistox(c) = parse_callisto_logx(c, callisto(c)...)

pzip(a) = Point.(a)
tzip(a::Series) = pzip(tuples(a))
tzip(a::Array)  = pzip.(tuples.(a))

plot(p, f; kwargs...) = save(drawplot(tzip(p) ; kwargs...), f)
save(d, fn) =  pk.save(d, plotpath(fn), 4)

getopts(; kw...) = kw

function setup()
    Random.seed!(1)
    num_nodes = 3
    typical_freq = 0.2 # GHz
    max_ppm = 100
    freqs = typical_freq * (1 .+  rand(-max_ppm:max_ppm, num_nodes) /  1e6)
    errors =  [Error(a) for a in freqs]
    
    return getopts(; topology = ["triangle"],
                   latency = 200,
                   errors,
                   poll_period = 1e5 * typical_freq,
                   control_delay = 10,
                   tmax = 2e8,
                   base_freq = typical_freq,
                   kp = 2e-8,
                   ki = 1e-15)
end


function pump_controller(; freq_step = 5e-7, poll_period = 1, base_freq = 1, kp =1, ki = 1, kw...)
    controller_init = () -> (0.0, 0.0)
    function controller_next(i, (xi, z), ml)
        r = sum(a[2] for a in ml)
        next_xi =  xi + poll_period*r/base_freq
        correction_des =  kp * r +  ki * next_xi
        pump = sign(correction_des - freq_step * z)
        next_z = z + pump
        correction = freq_step * next_z
        return (next_xi, next_z), correction
    end
    return controller_init, controller_next
end



function main()
    opts1 = setup()

    # first run without a pulsing controller
    c = CalOpts(; opts1...)
    results = callistox(c)

    # now run with the pulsing controller
    controller_init, controller_next = pump_controller(; opts1...)
    cp = CalOpts(; opts1..., controller_init, controller_next)
    resultsp = callistox(cp)

    # plots
    plot(results.freq,  "pu_simple_freqs1.pdf")
    plot(results.freq,  "pu_simple_freqs2.pdf"; xmin = 10e6, xmax=21e6)
    plot(results.mocc,  "pu_simple_mocc1.pdf")

    plot(resultsp.freq,  "pu_pulsed_freqs1.pdf")
    plot(resultsp.freq,  "pu_pulsed_freqs2.pdf"; xmin = 10e6, xmax=21e6)
    plot(resultsp.mocc,  "pu_pulsed_mocc1.pdf")

end


end

