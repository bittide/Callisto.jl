
module LogData

export Log, llog, llognext, num_allocated_records, get_records, get_fields, make_tuples

# row always points to an empty slot
mutable struct Log
    nfields
    num_preallocated_records
    current_record
    data
    header
end

function Log(nfields, header)
    #num_preallocated_records = 1000000
    num_preallocated_records = 1000
    data = zeros(Float64, nfields, num_preallocated_records)
    log = Log(nfields, num_preallocated_records, 1, data, header)
end

function llognext(logobj)
    logobj.current_record += 1
    if logobj.current_record > num_allocated_records(logobj)
        logobj.data = [logobj.data  zeros(logobj.nfields, logobj.num_preallocated_records)]
    end
end


function llog(logobj, ltype, a, b, c, d, e, f, g)
    record = logobj.current_record
    logobj[1, record] = ltype
    logobj[2, record] = a 
    logobj[3, record] = b 
    logobj[4, record] = c 
    logobj[5, record] = d 
    logobj[6, record] = e 
    logobj[7, record] = f
    logobj[8, record] = g
    llognext(logobj)
    return nothing
end

function llog(logobj, ltype, a, b, c, d, e, f, g, h)
    record = logobj.current_record
    logobj[1, record] = ltype
    logobj[2, record] = a 
    logobj[3, record] = b 
    logobj[4, record] = c 
    logobj[5, record] = d 
    logobj[6, record] = e 
    logobj[7, record] = f
    logobj[8, record] = g
    logobj[9, record] = h
    llognext(logobj)
    return nothing
end

function llog(logobj, ind, dat)
    record = logobj.current_record
    logobj[ind, record] = dat
    return nothing
end


Base.getindex(a::Log, field, record) = a.data[field, record]
Base.setindex!(a::Log, v, field, record) = (a.data[field, record] = v)
num_allocated_records(a::Log) = size(a.data, 2)
namedtuple(names, values) = NamedTuple{tuple(Symbol.(names)...)}(values)

function make_tuples(logobj)
    U = logobj.data[:, 1:logobj.current_record-1]
    U3 = [namedtuple(logobj.header,x) for x in eachcol(U)]
    return U3
end

function get_records(U3, f)
    #U2 = hcat([x  for x in eachcol(U) if f(x)]...)
    D2 = [x for x in U3 if f(x)]
end


function get_fields(logobj, fieldname)
    i = findfirst(x -> x == fieldname, logobj.header)
    x = logobj.data[i, 1:logobj.current_record-1]
end

end
