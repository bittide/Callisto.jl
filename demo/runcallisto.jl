
module Run

##############################################################################
using Callisto
##############################################################################

using PlotKit
const pk = PlotKit

callistox(c) = fast_post_process(c, callisto(c)...)


pzip(a) = Point.(a)
tzip(a::Series) = pzip(tuples(a))
tzip(a::Array)  = pzip.(tuples.(a))


qsave(d, fn) =  pk.save(d, fn, 4)
plot(p, fname; kwargs...) =  qsave(draw(Chart(tzip(p); kwargs...)), fname)


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
