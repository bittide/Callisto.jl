
module Post


export parse_callisto_log, parse_callisto_logx, focused_callisto_info, get_freq
export getlatency, getlambda, getbeta, getomega, getdata
export CalPost, FastPost, make_frequencies, make_measurements, make_measured_occupancy, fast_post_process

import ..SimCore:    beta, gamma
import ..LogData:    make_tuples, get_records, make_records, get_field, make_fieldmap, num_records
import ..Piecewise:  PiecewiseConstant, Samples, delay, after, before, PiecewiseLinear


Base.@kwdef mutable struct CalPost
    c
    simlog
    tuples
    start_ind_by_node
    theta
    freq = nothing
    meas = nothing
    xi = nothing
    mocc = nothing
end


function make_fieldmap(cpost::CalPost)
    cpost.fieldmap = make_fieldmap(cpost.simlog)
end

function CalPost(c, simlog, theta)
    tuples = make_tuples(simlog)
    sort!(tuples, by = x -> x.nid)
    start_ind_by_node = [1]
    nid = 2
    for k=1:length(tuples)
        if tuples[k].nid == nid
            push!(start_ind_by_node, k)
            nid += 1
        end
    end
    push!(start_ind_by_node, length(tuples)+1)
    cpost = CalPost(; c, simlog, tuples, start_ind_by_node, theta)
end
        
get_data_for_node(cpost::CalPost, i) = cpost.tuples[cpost.start_ind_by_node[i]:cpost.start_ind_by_node[i+1]-1]


function make_frequencies(cpost::CalPost)
    n = cpost.c.graph.n
    cpost.freq = PiecewiseConstant[]
    for i=1:n
        X = get_data_for_node(cpost, i)
        fq = [a.freq for a in X]
        sk = [a.sk for a in X]
        push!(cpost.freq, PiecewiseConstant(sk, fq[1:end-1]))
    end
end

function make_measurements(cpost::CalPost)
    n = cpost.c.graph.n
    cpost.meas = Samples[]
    for i=1:n
        X = get_data_for_node(cpost, i)
        mdata = [a.meas for a in X]
        tk = [a.tk for a in X]
        sam = Samples(tk[3:end], mdata[3:end])
        push!(cpost.meas, sam)
    end
end

function make_controller_state(cpost::CalPost)
    n = cpost.c.graph.n
    cpost.xi = PiecewiseConstant[]
    for i=1:n
        X = get_data_for_node(cpost, i)
        xdata = [a.xi for a in X]
        tk = [a.tk for a in X]
        cs = PiecewiseConstant(tk[3:end], xdata[4:end])
        push!(cpost.xi, cs)
    end
end

function make_measured_occupancy(cpost::CalPost)
    m = cpost.c.graph.m
    cpost.mocc = Array{Union{Missing,Samples}, 1}(missing, m)
    for link in cpost.c.links
        tk = cpost.meas[link.dst].x
        mo = [beta(link, t, cpost.theta) for t in tk]
        cpost.mocc[link.id] = Samples(tk, mo)
    end
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
            #mo[r] =  beta(dst, src, e, t, theta, c.ugn, c.latency, c.gears)
            mo[r] =  beta(c.links[e], t, theta)
        end
        occ[e] = Samples(times, mo)
    end
    return occ
end


Base.@kwdef mutable struct FastPost
    c
    simlog
    sorted_data
    start_ind_by_node
    theta
    freq = nothing
    meas = nothing
    xi = nothing
    mocc = nothing
end

invert_header(simlog, fieldname) = findfirst(x -> x == fieldname, simlog.header)

function FastPost(c, simlog, theta)
    nid_row = invert_header(simlog, "nid")
    n = num_records(simlog)
    p = sortperm(simlog.data[nid_row,1:n])
    sorted_data = simlog.data[:,p]
    start_ind_by_node = [1]
    nid = 2
    for k=1:n
        if sorted_data[nid_row, k] == nid
            push!(start_ind_by_node, k)
            nid += 1
        end
    end
    push!(start_ind_by_node, n+1)
    fpost = FastPost(; c, simlog, sorted_data, start_ind_by_node, theta)
end

function get_all_data(fpost::FastPost, fieldname, nid)
    field_row = invert_header(fpost.simlog, fieldname)
    return fpost.sorted_data[field_row, fpost.start_ind_by_node[nid]:fpost.start_ind_by_node[nid+1]-1]
end


function make_frequencies(fpost::FastPost)
    n = fpost.c.graph.n
    fpost.freq = PiecewiseConstant[]
    for i=1:n
        fq = get_all_data(fpost, "freq", i)
        sk = get_all_data(fpost, "sk", i)
        push!(fpost.freq, PiecewiseConstant(sk, fq[1:end-1]))
    end
end

    
function make_measurements(fpost::FastPost)
    n = fpost.c.graph.n
    fpost.meas = Samples[]
    for i=1:n
        mdata = get_all_data(fpost, "meas", i)
        tk = get_all_data(fpost, "tk", i)
        sam = Samples(tk[3:end], mdata[3:end])
        push!(fpost.meas, sam)
    end
end

function make_measured_occupancy(fpost::FastPost)
    m = fpost.c.graph.m
    fpost.mocc = Array{Union{Missing,Samples}, 1}(missing, m)
    for link in fpost.c.links
        tk = fpost.meas[link.dst].x
        mo = [beta(link, t, fpost.theta) for t in tk]
        fpost.mocc[link.id] = Samples(tk, mo)
    end
end
    
function fast_post_process(c, simlog, theta)
    fpost = FastPost(c, simlog, theta)
    make_frequencies(fpost)
    make_measurements(fpost)
    make_measured_occupancy(fpost)
    return fpost
end

    
function parse_callisto_logx(c, simlog, theta)
    cpost = CalPost(c, simlog, theta)
    make_frequencies(cpost)
    make_measurements(cpost)
    make_controller_state(cpost)
    make_measured_occupancy(cpost)
    return cpost
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



getlatency(c) = [l.latency for l in c.links]
getlambda(c) = [l.ugn for l in c.links]
getomega(t, xc) = [f(t) for f in xc.freq]
getbeta(t, c, xc) = [beta(l, t, xc.theta) for l in c.links]
getgamma(t, c, xc) = [gamma(l, t, xc.theta) for l in c.links]

function getbsdz(c)
    B = c.graph.incidence
    S = Int.(B .> 0)
    D = Int.(B .< 0)
    Z = c.graph.fundamental_cycles
    return B, S, D, Z
end

function getdata(t, c, xc)
    B, S, D, Z = getbsdz(c)
    d = (c = c, xc = xc, t = t,
         B = B, S = S, D = D, Z = Z,
         lambda = getlambda(c),
         latency = getlatency(c),
         beta0  = getbeta(0, c, xc),
         gamma0 = getgamma(0, c, xc),
         omega0 = getomega(0, xc),
         beta   = getbeta(t, c, xc),
         gamma  = getgamma(t, c, xc),
         omega  = getomega(t, xc))
    return d
end




end

