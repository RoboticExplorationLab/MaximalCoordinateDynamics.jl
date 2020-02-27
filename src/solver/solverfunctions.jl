@inline function setDandΔs!(diagonal::DiagonalEntry, body::Body, mechanism::Mechanism)
    diagonal.D = ∂dyn∂vel(body, mechanism.dt)
    diagonal.Δs = dynamics(body, mechanism)
    return
end

@inline function setDandΔs!(diagonal::DiagonalEntry{T,N}, eqc::EqualityConstraint, mechanism::Mechanism) where {T,N}
    diagonal.D = @SMatrix zeros(T, N, N)
    # μ = 1e-5
    # d.D = SMatrix{N,N,T,N*N}(μ*I) # TODO Positiv because of weird system? fix generally
    diagonal.Δs = g(eqc, mechanism)
    return
end

@inline function setLU!(offdiagonal::OffDiagonalEntry, bodyid::Int64, eqc::EqualityConstraint, mechanism)
    offdiagonal.L = -∂g∂pos(eqc, bodyid, mechanism)'
    offdiagonal.U = ∂g∂vel(eqc, bodyid, mechanism)
    return
end

@inline function setLU!(offdiagonal::OffDiagonalEntry, eqc::EqualityConstraint, bodyid::Int64, mechanism)
    offdiagonal.L = ∂g∂vel(eqc, bodyid, mechanism)
    offdiagonal.U = -∂g∂pos(eqc, bodyid, mechanism)'
    return
end

@inline function setLU!(offdiagonal::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    offdiagonal.L = @SMatrix zeros(T, N2, N1)
    offdiagonal.U = offdiagonal.L'
    return
end

@inline function updateLU!(offdiagonal::OffDiagonalEntry, diagonal::DiagonalEntry)
    Dinv = diagonal.Dinv
    offdiagonal.L = offdiagonal.L * Dinv
    offdiagonal.U = Dinv * offdiagonal.U
    return
end

@inline function updateD!(diagonal::DiagonalEntry, childdiagonal::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.D -= fillin.L * childdiagonal.D * fillin.U
    return
end

function invertD!(diagonal::DiagonalEntry)
    diagonal.Dinv = inv(diagonal.D)
    return
end

@inline function LSol!(diagonal::DiagonalEntry, child::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.Δs -= fillin.L * child.Δs
    return
end

function DSol!(diagonal::DiagonalEntry)
    diagonal.Δs = diagonal.Dinv * diagonal.Δs
    return
end

@inline function USol!(diagonal::DiagonalEntry, parent::DiagonalEntry, fillin::OffDiagonalEntry)
    diagonal.Δs -= fillin.U * parent.Δs
    return
end


function factor!(graph::Graph, ldu::SparseLDU)
    for id in graph.dfslist
        childs = directchildren(graph, id)
        for childid in childs
            offdiagonal = getentry(ldu, (id, childid))
            updateLU!(offdiagonal, getentry(ldu, childid))
        end

        diagonal = getentry(ldu, id)

        for childid in directchildren(graph, id)
            updateD!(diagonal, getentry(ldu, childid), getentry(ldu, (id, childid)))
        end
        invertD!(diagonal)
    end
end

function solve!(mechanism)
    ldu = mechanism.ldu
    graph = mechanism.graph
    dfslist = graph.dfslist

    for id in dfslist
        diagonal = getentry(ldu, id)

        for childid in directchildren(graph, id)
            LSol!(diagonal, getentry(ldu, childid), getentry(ldu, (id, childid)))
        end
    end

    for id in graph.rdfslist
        diagonal = getentry(ldu, id)

        DSol!(diagonal)

        for parentid in predecessors(graph, id)
            USol!(diagonal, getentry(ldu, parentid), getentry(ldu, (parentid, id)))
        end
    end
end

@inline function s0tos1!(component::Component)
    component.s1 = component.s0
    return
end

@inline function s1tos0!(component::Component)
    component.s0 = component.s1
    return
end


@inline function normΔs(component::Component)
    d = component.s1 - component.s0
    return dot(d, d)
end
