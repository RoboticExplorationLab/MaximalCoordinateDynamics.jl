abstract type Entry{T} end

mutable struct DiagonalEntry{T,N,N²} <: Entry{T}
    D::SMatrix{N,N,T,N²}
    Dinv::SMatrix{N,N,T,N²}
    Δs::SVector{N,T}

    function DiagonalEntry{T,N}() where {T,N}
        N² = N^2
        D = @SMatrix zeros(T, N, N)
        Dinv = @SMatrix zeros(T, N, N)
        Δs = @SVector zeros(T, N)

        new{T,N,N²}(D, Dinv, Δs)
    end
end

mutable struct OffDiagonalEntry{T,N1,N2,N1N2} <: Entry{T}
    L::SMatrix{N2,N1,T,N1N2}
    U::SMatrix{N1,N2,T,N1N2}

    function OffDiagonalEntry{T,N1,N2}() where {T,N1,N2}
        L = @SMatrix zeros(T, N2, N1)
        U = @SMatrix zeros(T, N1, N2)
        
        new{T,N1,N2,N1 * N2}(L, U)
    end
end


struct SparseLDU{T}
    diagonals::UnitDict{Base.OneTo{Int64},DiagonalEntry{T}}
    offdiagonals::Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}

    function SparseLDU(graph::Graph{N},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}},
        bdict::Dict,eqdict::Dict) where {T,N}

        diagonals = DiagonalEntry{T}[]
        for body in bodies
            push!(diagonals, DiagonalEntry{T,length(body)}())
        end
        for constraint in eqcs
            push!(diagonals, DiagonalEntry{T,length(constraint)}())
        end
        diagonals = UnitDict(diagonals)

        offdiagonals = Dict{Tuple{Int64,Int64},OffDiagonalEntry{T}}()
        for id in graph.dfslist
            haskey(bdict, id) ? node = bodies[bdict[id]] : node = eqcs[eqdict[id]]
            N1 = length(node)

            for cid in directchildren(graph, id)
                haskey(bdict, cid) ? cnode = bodies[bdict[cid]] : cnode = eqcs[eqdict[cid]]
                N2 = length(cnode)

                offdiagonals[(id, cid)] = OffDiagonalEntry{T,N2,N1}()
            end
        end

        new{T}(diagonals, offdiagonals)
    end
end

@inline getentry(ldu::SparseLDU, id::Int64) = ldu.diagonals[id]
@inline getentry(ldu::SparseLDU, ids::Tuple{Int64,Int64}) = ldu.offdiagonals[ids]
