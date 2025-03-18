
module FastPostProc

import ..SimCore:    beta
import ..LogData:    num_records
import ..Piecewise:  PiecewiseConstant, Samples

export  FastPost, fast_post_process, quick_post


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


invert_header(simlog, fieldname) = findfirst(x -> x == fieldname, simlog.header)



function get_all_data(fpost::FastPost, fieldname, nid; startoffset=0, endoffset=0)
    field_row = invert_header(fpost.simlog, fieldname)
    return fpost.sorted_data[field_row, (startoffset + fpost.start_ind_by_node[nid]):fpost.start_ind_by_node[nid+1]-1+endoffset]
end


function make_frequencies(fpost::FastPost)
    n = fpost.c.graph.n
    fpost.freq = PiecewiseConstant[]
    for i=1:n
        fq = get_all_data(fpost, "freq", i; endoffset=-1)
        sk = get_all_data(fpost, "sk", i)
        push!(fpost.freq, PiecewiseConstant(sk, fq))
    end
end


    
function make_controller_state(fpost::FastPost)
    n = fpost.c.graph.n
    fpost.xi = PiecewiseConstant[]
    for i=1:n
        xdata = get_all_data(fpost, "xi", i)
        tk = get_all_data(fpost, "tk", i)
        sam = PiecewiseConstant(tk[3:end], xdata[4:end])
        push!(fpost.xi, sam)
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
        #tk = fpost.meas[link.dst].x
        tk = get_all_data(fpost, "tk", link.dst; startoffset=2)
        mo = [beta(link, t, fpost.theta) for t in tk]
        fpost.mocc[link.id] = Samples(tk, mo)
    end
end


function fast_post_process(c, simlog, theta)
    fpost = FastPost(c, simlog, theta)
    make_frequencies(fpost)
    #make_measurements(fpost)
    make_measured_occupancy(fpost)
    return fpost
end


##############################################################################
# For those who are extra hasty, use this.
# This is not approximate, but it only provides equidistant sampled values.


#
# Sample the value of the occupancy at a given list of times
#
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

#
# compute the sampled frequency by differentiating theta
#
function get_freq_samples(theta, tvals)
    freqs = []
    for t in tvals
        i = searchsortedlast(theta.x, t)
        if i == length(theta.x)
            i = i - 1
        end
        push!(freqs, (theta.y[i+1] - theta.y[i])/(theta.x[i+1] - theta.x[i]))
    end
    return freqs
end


#
# Sample the occupancy at N equidistant points between 0 and the
# max simulation time. Also sample the frequency, by differentiating
# theta, instead of using the logged value.
#
function quick_post(c, theta; N = 500)
    tmax = minimum(a.x[end] for a in theta)
    n = length(theta)
    tvals = collect(range(0, tmax, length = N))
    freq = Samples[]
    for i = 1:n
        fr = get_freq_samples(theta[i], tvals)
        push!(freq, Samples(tvals, fr))
    end
    occ = get_occ_samples(c, theta, tvals)
    qc = (freq = freq, occ = occ)
    return qc
end


end
