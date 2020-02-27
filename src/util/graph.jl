struct Graph{N}
    directchildren::Vector{Vector{Int64}}
    predecessors::Vector{Vector{Int64}}
    connections::Vector{Vector{Int64}}

    dfslist::SVector{N,Int64}
    rdfslist::SVector{N,Int64}

    dict::UnitDict{Base.OneTo{Int64},Int64}
    rdict::UnitDict{Base.OneTo{Int64},Int64}

    function Graph(origin::Origin,bodies::Vector{<:Body},
        eqconstraints::Vector{<:EqualityConstraint})

        oid = origin.id
        adjacency, dict = adjacencyMatrix(eqconstraints, bodies)
        dfsgraph, dfslist, loops = dfs(adjacency, dict, oid)

        adjacency = deleteat(adjacency, dict[oid])
        dfsgraph = deleteat(dfsgraph, dict[oid])
        dfslist = StaticArrays.deleteat(dfslist, length(dfslist))

        for (id, ind) in dict
            ind > dict[oid] && (dict[id] = ind - 1)
        end
        pop!(dict, oid)
        rdict = Dict(ind => id for (id, ind) in dict)

        for constraint in eqconstraints
            constraint.pid == oid && (constraint.pid = nothing)
        end

        N = length(dict)

        adjacency = convert(Vector{SVector{N,Bool}}, adjacency)
        dfsgraph = convert(Vector{SVector{N,Bool}}, dfsgraph)

        dirs = directchildren(dfslist, dfsgraph, dict)
        preds = predecessors(dfslist, dfsgraph, dict)
        cons = connections(dfslist, adjacency, dict)

        dict = UnitDict(dict)
        rdict = UnitDict(rdict)

        new{N}(dirs, preds, cons, dfslist, reverse(dfslist), dict, rdict)
    end
end

function adjacencyMatrix(eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body})
    A = zeros(Bool, 0, 0)
    dict = Dict{Int64,Int64}()
    n = 0

    for constraint in eqconstraints
        cid = constraint.id
        A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
        dict[cid] = n += 1
        for bodyid in unique([constraint.pid;constraint.bodyids])
            if !haskey(dict, bodyid)
                A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
                dict[bodyid] = n += 1
            end
            A[dict[cid],dict[bodyid]] = true
        end
    end
    for body in bodies # add unconnected bodies
        if !haskey(dict, body.id)
            A = [A zeros(Bool, n, 1); zeros(Bool, 1, n) zero(Bool)]
            dict[body.id] = n += 1
        end
    end

    A = A .| A'
    convert(Matrix{Bool}, A), dict
end

function dfs(adjacency::Matrix, dict::Dict, originid::Int64)
    N = size(adjacency)[1]
    dfsgraph = zeros(Bool, N, N)
    dfslist = zeros(Int64, N)
    visited = zeros(Bool, N)
    loops = Vector{Int64}[]
    index = N

    dfslist[index] = originid
    visited[dict[originid]] = true
    dfs!(adjacency, dfsgraph, dict, dfslist, visited, loops, N, originid, -1)
    loops = loops[sortperm(sort.(loops))[1:2:length(loops)]] # removes double entries of loop connections and keeps the first found pair

    return dfsgraph, convert(SVector{N}, dfslist), loops
end

function dfs!(A::Matrix, Adfs::Matrix, dict::Dict, list::Vector, visited::Vector, loops::Vector{Vector{Int64}}, index::Int64, currentid::Int64, parentid::Int64)
    i = dict[currentid]
    for (childid, j) in dict
        if A[i,j] && parentid != childid # connection from i to j in adjacency && not a direct connection back to the parent
            if visited[j]
                push!(loops, [childid,currentid]) # childid is actually a predecessor of currentid since it's a loop
            else
                index -= 1
                list[index] = childid
                visited[j] = true
                Adfs[i,j] = true
                index = dfs!(A, Adfs, dict, list, visited, loops, index, childid, currentid)
            end
        end
    end
    return index
end

function parent(dfsgraph::Matrix, dict::Dict, childid::Int64) where {N,T}
    j = dict[childid]
    for (parentid, i) in dict
        dfsgraph[i,j] && (return parentid)
    end
    return -1
end

# this is done in order!
function directchildren(dfslist, dfsgraph, dict::Dict)
    N = length(dfslist)
    dirs = [Int64[] for i = 1:N]
    for i = 1:N
        for cid in dfslist
            dfsgraph[i][dict[cid]] && push!(dirs[i], cid)
        end
    end

    return dirs
end

# this is done in reverse order (but this is not really important for predecessors)
function predecessors(dfslist, pattern, dict::Dict)
    N = length(dfslist)
    preds = [Int64[] for i = 1:N]
    for i = 1:N
        for cid in reverse(dfslist)
            pattern[dict[cid]][i] && push!(preds[i], cid)
        end
    end

    return preds
end

# this is done in order (but this is not really important for connections)
function connections(dfslist, adjacency, dict::Dict)
    N = length(dfslist)
    cons = [Int64[] for i = 1:N]
    for i = 1:N
        for cid in dfslist
            adjacency[i][dict[cid]] && push!(cons[i], cid)
        end
    end

    return cons
end

@inline directchildren(graph, id::Int64) = graph.directchildren[graph.dict[id]]
@inline predecessors(graph, id::Int64) = graph.predecessors[graph.dict[id]]
@inline connections(graph, id::Int64) = graph.connections[graph.dict[id]]
