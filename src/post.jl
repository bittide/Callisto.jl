
module Post


export focused_callisto_info

import Piecewise: after, before, delay


##############################################################################
# For those who want all the details, use this

# buffer at node i associated with link to node j
# find occupancy as a piecewise constant function
function beta_piecewise(c, link, theta, tmin, tmax)
    thj1 = delay(theta[link.src], link.latency)
    thj = before(after(thj1, tmin), tmax)
    thi = before(after(theta[link.dst], tmin), tmax)
    p = floor(link.gear*thj) - floor(link.gear*thi)
    return p + link.ugn
end


# things that take a long time we only apply to a time range
function get_occ(c, theta, tmin, tmax)
    occ = [beta_piecewise(c, link, theta, tmin, tmax) for link in c.model.links]
    return occ
end

function focused_callisto_info(c, xc, tmin, tmax)
    occ = get_occ(c, xc.theta, tmin, tmax)
    fc = (occ = occ, )
    return fc
end



end

