
module Run

##############################################################################
using Callisto
##############################################################################

using PlotKit
const pk = PlotKit

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)
callistox(c) = parse_callisto_logx(c, callisto(c)...)


pzip(a) = Point.(a)
tzip(a::Series) = pzip(tuples(a))
tzip(a::Array)  = pzip.(tuples.(a))

plot(p, f; kwargs...) = save(drawplot(tzip(p) ; kwargs...), f)
save(d, fn) =  pk.save(d, plotpath(fn), 4)



# simplest run
function main()
    c = CalOpts(; kp = 2e-8)
    xc = callistox(c)
    plot(xc.freq,  "freqs.pdf")
    plot(xc.theta, "phase.png")
    plot(xc.mocc,  "mocc.pdf"; xmax=500e6)

    fc = focused_callisto_info(c, xc, 125e6, 125e6 + 1e2)
    plot(fc.occ, "occupancy.pdf")

    return xc
end


end
