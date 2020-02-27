# Newton with line search 
function newton!(mechanism::Mechanism{T,Nl,0}; ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false) where {T,Nl}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(mechanism) # x̂1 for each body and constraint

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        normΔs1 = normΔs(mechanism)

        foreach(s1tos0!, bodies)
        foreach(s1tos0!, eqcs)

        if normΔs1 < ε && normf1 < ε
            # warning && (@info string("Newton iterations: ",n))
            return
        else
            normf0 = normf1
        end
    end

    warning && (@info string("newton_ip! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
    return
end

function lineSearch!(mechanism::Mechanism{T,N,0}, normf0;iter = 10, warning::Bool = false) where {T,N}
    normf1 = normf0
    scale = 0
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for n = Base.OneTo(iter + 1)
        for body in bodies
            lineStep!(body, getentry(ldu, body.id), scale)
        end
        for eqc in eqcs
            lineStep!(eqc, getentry(ldu, eqc.id), scale)
        end

        normf1 = normf(mechanism)
        if normf1 >= normf0
            scale += 1
        else
            return normf1
        end
    end

    warning && (@info string("lineSearch! did not converge. n = ", iter, "."))
    return normf1
end


@inline function lineStep!(component::Component, diagonal::DiagonalEntry, scale)
    component.s1 = component.s0 - 1 / (2^scale) * diagonal.Δs
    return
end