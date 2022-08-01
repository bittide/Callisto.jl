
module Post


export parse_callisto_log, parse_callisto_logx, focused_callisto_info, get_freq

import ..SimCore:    beta
import ..LogData:    make_tuples, get_records
import ..Piecewise:  PiecewiseConstant, Samples, delay, after, before, PiecewiseLinear


function get_freq(c, Tdata)
    n = c.graph.n
    freq = PiecewiseConstant[]
    for i=1:n
        X = get_records(Tdata, x -> x.nid == i)
        fq = [a.freq for a in X]
        sk = [a.sk for a in X]
        push!(freq, PiecewiseConstant(sk, fq[1:end-1]))
    end
    return freq
end

function get_occ_samples(c, theta, times)   
    # construct occupancies at measurement times
    n = c.graph.n
    m = c.graph.m
    occ = Array{Union{Missing,Samples}, 1}(missing, m)
    for e = 1:m
        src, dst = c.graph.edges[e]
        mo = zeros(length(times))
        for r = 1:length(times)
            t = times[r]
            mo[r] =  beta(dst, src, e, t, theta, c.ugn, c.latency, c.gears)
        end
        occ[e] = Samples(times, mo)
    end
    return occ
end

function parse_callisto_log(c, Tdata, theta)
    p = c.poll_period
    d = c.control_delay
    n = c.graph.n
    m = c.graph.m
    
    meas = Samples[]
    for i=1:n
        X = get_records(Tdata, x -> x.nid == i)
        mdata = [a.meas for a in X]
        tk = [a.tk for a in X]
        sam = Samples(tk[3:end], mdata[3:end])
        push!(meas, sam)
    end

    # at time tk the controller state
    # changes from xi(k) to xi(k+1)
    xi = PiecewiseConstant[]
    for i=1:n
        X = get_records(Tdata, x -> x.nid == i)
        xdata = [a.xi for a in X]
        tk = [a.tk for a in X]
        #cs = PiecewiseConstant(tk[3:end], xdata[3:end-1])
        cs = PiecewiseConstant(tk[3:end], xdata[4:end])
        push!(xi, cs)
    end

    
    # construct occupancies at measurement times
    mocc = Array{Union{Missing,Samples}, 1}(missing, m)
    for link in c.links
        tk = meas[link.dst].x
        mo = [beta(link, t, theta) for t in tk]
        mocc[link.id] = Samples(tk, mo)
    end
    return meas, xi, mocc
end


function parse_callisto_logx(c, simlog, theta)
    Tdata = make_tuples(simlog)
    freq = get_freq(c, Tdata)
    meas, xi, mocc = parse_callisto_log(c, Tdata, theta)  
    adjusted_freqs = [ 10^9*(f + (-1)) for f in freq]

    return (simlog=simlog, theta=theta, freq=freq, meas=meas,
            afreq=adjusted_freqs, xi=xi, mocc=mocc)
end


# things that take a long time we only apply to a time range
function get_occ(c, theta, tmin, tmax)
    occ = [beta_piecewise(c, link, theta, tmin, tmax) for link in c.links]
    return occ
end

function get_delta(c, theta)
    delta = Array{Union{Missing, PiecewiseLinear}, 2}(missing, n, n)
    for i=1:n
        for j in c.graph.adjacent_nodes[i]
            delta[i,j] = delay(theta[j], c.latency[i,j]) - theta[i] + c.ugn[i,j]
        end
    end
    return delta
end


# buffer at node i associated with link to node j
# find occupancy as a piecewise constant function
function beta_piecewise(c, link, theta, tmin, tmax)
    thj1 = delay(theta[link.src], link.latency)
    thj = before(after(thj1, tmin), tmax)
    thi = before(after(theta[link.dst], tmin), tmax)
    p = floor(link.gear*thj) - floor(link.gear*thi)
    return p + link.ugn
end

function get_adjusted_theta(c, theta,  tmin, tmax)
    n = c.graph.n
    proj(x) = before(after(x, tmin),tmax)
    adj_freq = minimum( (theta[i](tmax) - theta[i](tmin)) / (tmax-tmin) for i = 1:n)
    adj_freq = adj_freq - 0.1/(tmax-tmin)
    epoch = -3*maximum([a.latency for a in c.links])
    adj_theta = PiecewiseLinear[]
    for i=1:n
        rateadj = PiecewiseLinear([epoch, tmax],[epoch*adj_freq, tmax*adj_freq])
        ath = proj(theta[i]) - rateadj
        ath = ath - floor(theta[i](tmin))
        push!(adj_theta, ath)
    end
    return adj_theta
end

function focused_callisto_info(c, xc, tmin, tmax)
    occ = get_occ(c, xc.theta, tmin, tmax)
    adj_theta = get_adjusted_theta(c, xc.theta, tmin, tmax)
    fc = (occ = occ, adj_theta = adj_theta)
    return fc
end




end

