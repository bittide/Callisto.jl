module Topology

using LinearAlgebra

export Graph, resistance, index_by_node


"""
    mutable struct Graph

A directed graph.

Fields are:
 - edges            list of edges. Each edge is a Tuple of the form `(i,j)`
                    where `i` and `j` are integers in 1,...,n.
 - n                number of nodes
 - m                number of edges
 - incidence        n by m incidence matrix
 - incoming_nodes   incoming_nodes[i] is a list of incoming neighbors of i
 - adjacent_edges   adjacent_edges[i] is a list of edges that connect to node i
 - topology         tuple

The source data here consists of `edges` and `n`. Apart from `topology`,
everything else is computed from those two quantities.
"""
mutable struct Graph
    edges               # list of named tuples, numbered 1,...,m, where m = length(edges)
    n                   # number of nodes. Nodes are numbered 1,...,n
    m                   # number of edges
    incidence           # n by m incidence matrix
    fundamental_cycles  # m by (m-n+1) matrix whose cols are the fcycles
    incoming_nodes      # incoming_nodes[i] is a list of incoming neighbors of i
    adjacent_edges      # adjacent_edges[i] is a list of edges that connect to node i
    outgoing_edges      # outgoing_edges[i] is a list of edges that depart from node i
    incoming_edges      # incoming_edges[i] is a list of edges which arrive at node i
    topology            # tuple
    layout              # a list of tuples giving coordinates in the plane
end


#
# TODO: make edges be named tuples, with keys src, dst.
#

"""
    Graph(topology)

Returns a graph with the specified topology.

The `topology` argument is a Tuple, for example `("mesh", 3, 4)`.
The graphs are undirected but oriented, so if `(1,2)` is
an edge then `(2,1)` will not be.
The first element is a string, giving the name of the class
of topologies, and the remaining elements are parameters.
Possible topologies include

 - ("triangle")          The simplest graph
 - ("diamond")           The second simplest
 - ("mesh", 3, 4)        A mesh graph with 3x4 nodes
 - ("full", n)           Complete graph
 - ("line", n)           Line graph with n nodes
 - ("hypercube", d)      Hypercube graph of order d
 - ("star", d)           Star graph with n+1 vertices
 - ("torus2d", x, y)     Torus 2d with side lengths x,y
 - ("torus3d", x, y, z)  Torus 3d with side lengths x,y,z
 - ("tree", d, c)        Tree with depth d and c children per level

"""
function Graph(topology; bidirectional = false, treeedgesfirst = false)
    edges, n, layout = get_topology(topology)
    if bidirectional
        edges = make_bidirectional(edges)
    end
    return Graph(edges, n; layout, topology, treeedgesfirst)
end


"""
    Graph(edges, n)

Constructor. `edges` is a list of Tuples of node ids, `n` is the number of nodes.
"""
function Graph(edges, n; layout = nothing, topology = ("unknown"), treeedgesfirst = false)
    m = length(edges)
    if treeedgesfirst
        edges = edges[find_spanning_tree(edges, n)]
    end
    incoming_nodes = get_incoming_nodes(edges, n)
    adjacent_edges = get_adjacent_edges(edges, n)
    incoming_edges = get_incoming_edges(edges, n)
    outgoing_edges = get_outgoing_edges(edges, n)
    incidence = get_incidence(edges, n)
    if treeedgesfirst
        fundamental_cycles = find_fundamental_cycles(incidence)
    else
        fundamental_cycles = nothing
    end
    g = Graph(edges, n, m, incidence, fundamental_cycles, incoming_nodes,
              adjacent_edges, outgoing_edges, incoming_edges, topology, layout)
end


"""
    find_fundamental_cycles(B)

Returns a matrix whose columns are incidence vectors
for the fundamental cycles of the graph. In particular
they form an integral basis for null(B). Note that B must
be sorted so that the first n-1 columns are the edges of
a spanning tree.
"""
function find_fundamental_cycles(B)
    n, m = size(B)
    B11 = B[1:n-1, 1:n-1]
    B12 = B[1:n-1, n:m]
    N = B11 \ B12
    Z = [-Int.(round.(N)); I(m-n+1)]
    return Z
end



"""
    find_spanning_tree(edges, n)

Return a permutation vector p of length m such that p[1:n-1] are
the edge ids of a spanning tree.

"""
function find_spanning_tree(edges, n)
    m = length(edges)
    treeedges = []
    nontreeedges = []
    treenodes = Set()
    for i = 1:m
        e = edges[i]
        if ! (e.src in treenodes) || !(e.dst in treenodes)
            push!(treeedges, i)
            push!(treenodes, e.src)
            push!(treenodes, e.dst)
        else
            push!(nontreeedges, i)
        end
    end
    p = vcat(treeedges, nontreeedges)
    return p
end

"""
    get_incoming_nodes(edges, n)

Return the neighbors of every node.

Return a list `incoming_nodes` such that `incoming_nodes[i]` is the list
of incoming neighbors of node `i`.
"""
function get_incoming_nodes(edges, n)
    incoming_nodes = [Int64[] for i=1:n]
    for (src, dst) in edges
        #push!(incoming_nodes[src], dst)
        push!(incoming_nodes[dst], src)
    end
    return incoming_nodes
end

"""
    get_incidence(edges, n)

Return the incidence matrix of the graph.

The incidence matrix `B` is `n` by `m`, where `n` is the number of nodes
and `m` is the number of edges. It is defined by

    B[i,j] = 1  if edge j starts at link i
    B[i,j] = -1 if edge j ends at link i
"""
function get_incidence(edges, n)
    m = length(edges)
    B = zeros(Int64, n, m)
    for (j, (src, dst)) in enumerate(edges)
        B[src, j] = 1
        B[dst, j] =-1
    end
    return B
end

"""
    get_adjacent_edges(edges, n)

Return the list of edges adjacent to each node

For each node `i`,  `adjacent_edges[i]` is a list of the edge ids
that either start or end at `i`
"""
function get_adjacent_edges(edges, n)
    adjacent_edges = [Int64[] for i=1:n]
    for (i, (src, dst)) in enumerate(edges)
        push!(adjacent_edges[src], i)
        push!(adjacent_edges[dst], i)
    end
    return adjacent_edges
end


"""
    get_outgoing_edges(edges, n)

Return the list of edges outgoing from each node

For each node `i`,  `outgoing_edges[i]` is a list of the edge ids
that either start or end at `i`
"""
function get_outgoing_edges(edges, n)
    outgoing_edges = [Int64[] for i=1:n]
    for (i, (src, dst)) in enumerate(edges)
        push!(outgoing_edges[src], i)
    end
    return outgoing_edges
end



"""
    get_incoming_edges(edges, n)

Return the list of edges incoming to each node

For each node `i`,  `incoming_edges[i]` is a list of the edge ids
that end at `i`
"""
function get_incoming_edges(edges, n)
    incoming_edges = [Int64[] for i=1:n]
    for (i, (src, dst)) in enumerate(edges)
        push!(incoming_edges[dst], i)
    end
    return incoming_edges
end




edge(s, d)  = (src = s, dst=d)


function meshlayout(nx, ny)
     n = nx * ny
     xlist = Tuple{Float64, Float64}[]
     for i = 1:n
         xc = (i-1) % nx
         yc = div(i-1, nx)
         push!(xlist, (xc, yc))
     end
     return xlist
end

dist(a,b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
function hcube(i, d)
    s = bitstring(i)[end-d+1:end]
    z = [parse(Int, s[k]) for k=1:d]
    return 2z .- 1
end

"""
    get_topology(topology)

Returns a list of directed edges, and the number of nodes.
See `Graph` for details of the topology argument.
"""
function get_topology(topology)
    edges = Any[]
    if topology[1] == "triangle"
        edges = [edge(1,2), edge(1,3), edge(2,3)] # topology
        n = 3
        layout = [(0.5, -1/sqrt(12)), (-0.5, -1/sqrt(12)), (0.0, 1/sqrt(3))]

    elseif topology[1] == "diamond"
        edges = [edge(1,2), edge(2,3), edge(3,4), edge(4,2), edge(1,3)]
        layout = meshlayout(2,2)
        n = 4

    elseif topology[1] == "full"
        n = topology[2]
        edges = Any[]
        for i=1:n
            for j=1:i-1
                # j is src
                push!(edges, edge(j,i))
            end
        end
        layout = [(cos(2*pi*i/n), sin(2*pi*i/n)) for i=0:n-1]
        r = dist(layout[1], layout[2])
        layout = [(a/r,b/r) for (a,b) in layout]

    elseif topology[1] == "line"
        n = topology[2]
        edges = Any[]
        for i=2:n
            push!(edges, edge(i-1,i))
        end
        layout = meshlayout(n,1)

    elseif topology[1] == "hypercube"
        degree = topology[2]
        n = 2^degree
        edges = Any[]
        for i=0:n-1
            for b=0:degree-1
                if i & (1 << b) == 0
                    j = i | (1 << b)
                    push!(edges, edge(i+1, j+1))
                end
            end
        end
        if degree == 1
            layout = meshlayout(2,1)
        elseif degree == 2
            layout = meshlayout(2,2)
        elseif degree == 3
            x0 = meshlayout(2,2)
            x1 = [(a-0.5, b-0.5) for (a,b) in x0]
            x2 = [(3a,3b) for (a,b) in x1]
            layout = vcat(x1, x2)
        elseif degree == 4
            points4d = [hcube(i,degree) for i = 0:2^degree-1]
            q = sqrt(2)
            A = [sqrt(2)  0       -1  -1
                 0        sqrt(2)  1  -1 ]
            x1 = [A*z for z in points4d]
            layout = [(a[1],a[2]) for a in x1]
        else
            layout = nothing
        end

    elseif topology[1] == "mesh"
        mx = topology[2]
        my = topology[3]
        n = mx*my
        edges = Any[]
        for y=0:my-1
            for x=0:mx-1
                nbrs = Any[]
                if x < mx-1 && y < my-1
                    nbrs = [(x+1, y), (x, y+1)]
                elseif x == mx-1 && y < my-1
                    nbrs = [(x, y+1)]
                elseif x < mx-1 && y == my-1
                    nbrs = [(x+1,y)]
                end
                for a in nbrs
                    x1, y1 = a
                    push!(edges, edge(x+y*mx+1, x1+y1*mx+1))
                end
            end
        end
        layout = meshlayout(mx, my)

    elseif topology[1] == "star"
        hosts = topology[2]
        n = hosts + 1
        edges = Any[]
        for y = 1:hosts
            push!(edges, edge(1, y+1))
        end
        layout = [(cos(2*pi*i/hosts), sin(2*pi*i/hosts)) for i=0:hosts-1]
        pushfirst!(layout, (0.0, 0.0))

    elseif topology[1] == "torus2d"
        mx = topology[2]
        my = topology[3]
        n = mx*my
        edges = Any[]
        for x=0:mx-1
            for y=0:my-1
                nbrs = [(mod(x+1, mx), y), (x, mod(y+1, my))]
                for a in nbrs
                    x1, y1 = a
                    push!(edges, edge(x+y*mx+1, x1+y1*mx+1))
                end
            end
        end
        layout = nothing

    elseif topology[1] == "torus3d"
        mx = topology[2]
        my = topology[3]
        mz = topology[4]
        n = mx*my*mz
        edges = Any[]
        for x=0:mx-1
            for y=0:my-1
                for z=0:mz-1
                    nbrs = [(mod(x+1, mx), y, z),
                            (x, mod(y+1, my), z),
                            (x, y, mod(z+1,mz))]
                    for a in nbrs
                        x1, y1, z1 = a
                        push!(edges, edge(x+y*mx+z*mx*my+1, x1+y1*mx+z1*mx*my+1))
                    end
                end
            end
        end
        layout = nothing

    elseif topology[1] == "tree"
        height = topology[2]
        children = topology[3]
        edges = Any[]
        prev_nodes = 0
        for level=1:height-1
            prev_nodes += children^(level-1)
            for c=0:children^(level)-1
                child_id = c + prev_nodes
                parent_id = Int(floor(c/children)) + prev_nodes - children^(level-1)
                push!(edges, edge(parent_id+1, child_id+1))
            end
        end
        n = (1 - children^(height))/(1-children)
        layout = nothing
    end

    # At the moment all of the above topologies satisfy this property.
    # But one could have a topology where node n is isolated.
    # This test will then fail.
    @assert n == maximum(max(s,d) for (s,d) in edges)
    return edges, n, layout
end

"""
    index_by_node(edges, x)

Construct a matrix `X` containing the elements of `x` indexed by node.

The vector `x` is indexed by edge number, that is `x[e]` is associated
with edge `e`.  If edge `e` connects nodes `i` to `j`, then `X[i,j] = x[e]`.
If there is no edge from `i` to `j`, then `X[i,j] = missing`.
"""
function index_by_node(edges, x)
    n, m = get_dims(edges)
    X = Matrix{Union{Missing, eltype(x)}}(missing, n, n)
    for i = 1:m
        X[edges[i]...] = x[i]
    end
    return X
end

"""
    make_bidirectional(edges)

Construct the bidirectional graph given a unidirectional graph.
"""
function make_bidirectional(edges)
    reversed_edges = [ edge(dst, src) for (src, dst) in edges]
    bidirectional_edges = vcat(copy(edges), copy(reversed_edges))
    return bidirectional_edges
end

"""
    resistance(B)

Return the resistance matrix of a graph. `B` is the incidence matrix.
"""
function resistance(B)
    n,m = size(B)
    L = B * B';
    Q = pinv(L)
    R = zeros(n,n)
    for i=1:n
        for j=1:n
            R[i,j] = Q[i,i] + Q[j,j] - 2*Q[i,j]
        end
    end
    return R
end

end
