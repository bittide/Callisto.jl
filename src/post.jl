
module Post


export parse_callisto_logx, focused_callisto_info
export CalPost

import ..SimCore:    beta, gamma
import ..LogData:    make_tuples, get_records, make_records, get_field, make_fieldmap
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
   
    
function parse_callisto_logx(c, simlog, theta)
    cpost = CalPost(c, simlog, theta)
    make_frequencies(cpost)
    make_measurements(cpost)
    make_controller_state(cpost)
    make_measured_occupancy(cpost)
    return cpost
end



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
    occ = [beta_piecewise(c, link, theta, tmin, tmax) for link in c.links]
    return occ
end

function focused_callisto_info(c, xc, tmin, tmax)
    occ = get_occ(c, xc.theta, tmin, tmax)
    fc = (occ = occ, )
    return fc
end



end

