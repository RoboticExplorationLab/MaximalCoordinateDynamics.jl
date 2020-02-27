mutable struct Mechanism{T,N,Ni}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T
    No::Int64 # order of integrator, currently only No=2 (1st order) implemented

    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    graph::Graph{N}

    ldu::SparseLDU{T}
    storage::Storage{T}


    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
        eqcs::Vector{<:EqualityConstraint{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T


        resetGlobalID()

        Nb = length(bodies)
        Ne = length(eqcs)
        Ni = 0
        N = Nb + Ne
        steps = Int(ceil(tend / dt))

        currentid = 1

        bdict = Dict{Int64,Int64}()
        for (ind, body) in enumerate(bodies)
            push!(body.x, [body.x[1] for i = 1:No - 1]...)
            push!(body.q, [body.q[1] for i = 1:No - 1]...)
            push!(body.F, [body.F[1] for i = 1:No - 1]...)
            push!(body.τ, [body.τ[1] for i = 1:No - 1]...)

            for eqc in eqcs
                eqc.pid == body.id && (eqc.pid = currentid)
                for (ind, bodyid) in enumerate(eqc.bodyids)
                    if bodyid == body.id
                        eqc.bodyids = setindex(eqc.bodyids, currentid, ind)
                        eqc.constraints[ind].cid = currentid
                    end
                end
            end

            body.id = currentid
            currentid += 1

            bdict[body.id] = ind
        end

        eqdict = Dict{Int64,Int64}()
        for (ind, eqc) in enumerate(eqcs)
            eqc.id = currentid
            currentid += 1

            eqdict[eqc.id] = ind
        end

        normf = 0
        normΔs = 0

        graph = Graph(origin, bodies, eqcs)
        ldu = SparseLDU(graph, bodies, eqcs, bdict, eqdict)

        storage = Storage{T}(steps, Nb, Ne)

        bodies = UnitDict(bodies)
        eqcs = UnitDict((eqcs[1].id):(eqcs[Ne].id), eqcs)

        new{T,N,Ni}(tend, Base.OneTo(steps), dt, g, No, origin, bodies, eqcs, normf, normΔs, graph, ldu, storage)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(OriginConnection(origin, body)))
        end
        Mechanism(origin, bodies, eqc, tend = tend, dt = dt, g = g, No = No)
    end    
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T,N,0}) where {T,N}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies and ", length(M.eqconstraints), " constraints")
end

function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id, body) in pairs(mechanism.bodies)
        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), id, geteqconstraint(mechanism, cid), mechanism)
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, body, mechanism)
    end

    for node in mechanism.eqconstraints
        id = node.id

        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), node, cid, mechanism)
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, node, mechanism)
    end
end

@inline getbody(mechanism::Mechanism, id::Int64) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin
@inline geteqconstraint(mechanism::Mechanism, id::Int64) = mechanism.eqconstraints[id]
@inline getineqconstraint(mechanism::Mechanism, id::Int64) = mechanism.ineqconstraints[id]


@inline function normf(body::Body{T}, mechanism::Mechanism) where T
    f = dynamics(body, mechanism)
    return dot(f, f)
end

@inline function normf(c::EqualityConstraint, mechanism::Mechanism)
    f = g(c, mechanism)
    return dot(f, f)
end

@inline function GtλTof!(body::Body, eqc::EqualityConstraint, mechanism)
    body.f -= ∂g∂pos(eqc, body.id, mechanism)' * eqc.s1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(body, mechanism)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)

    return sqrt(mechanism.normf)
end


@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    # Allocates otherwise
    mechanism.normΔs += mapreduce(normΔs, +, mechanism.bodies)
    foreach(addNormΔs!, mechanism.eqconstraints, mechanism)

    return sqrt(mechanism.normΔs)
end


@inline function addNormf!(eqc::EqualityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(eqc, mechanism)
    return
end

@inline function addNormΔs!(component::Component, mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end


function saveToTraj!(mechanism::Mechanism, t)
    No = mechanism.No
    for (ind, body) in enumerate(mechanism.bodies)
        mechanism.storage.x[ind][t] = body.x[No]
        mechanism.storage.q[ind][t] = body.q[No]
    end
    for (ind, constraint) in enumerate(mechanism.eqconstraints)
        mechanism.storage.λ[ind][t] = constraint.s1
    end
end

@inline function updatePos!(body::Body, dt)
    x2 = body.x[2]
    q2 = body.q[2]
    body.x[1] = x2
    body.x[2] = x2 + getvnew(body) * dt
    body.q[1] = q2
    body.q[2] = dt / 2 * (Lmat(q2) * ωbar(body, dt))
    return
end


function simulate!(mechanism::Mechanism;save::Bool = false,debug::Bool = false)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    dt = mechanism.dt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)

    for i = mechanism.steps
        newton!(mechanism, warning = debug)
        save && saveToTraj!(mechanism, i)
        foreach(updatePos!, bodies, dt)

        # debug && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return
end


function plotθ(mechanism::Mechanism{T}, id) where T
    n = length(mechanism.bodies)
    θ = zeros(T, n, length(mechanism.steps))
    for i = 1:n
        qs = mechanism.storage.q[i]
        for (t, q) in enumerate(qs)
            θ[i,t] = angleaxis(q)[1] * sign(angleaxis(q)[2][1])
        end
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), θ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), θ[ind,:])
    end
    return p
end

function plotλ(mechanism::Mechanism{T}, id) where T
    n = sum(length.(mechanism.eqconstraints))
    λ = zeros(T, n, length(mechanism.steps))
    startpos = 1
    endpos = 0
    for i = 1:length(mechanism.eqconstraints)
        endpos = startpos + length(mechanism.eqconstraints[i]) - 1

        λs = mechanism.storage.λ[i]
        for (t, val) in enumerate(λs)
            λ[startpos:endpos,t] = val
        end

        startpos = endpos + 1
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), λ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), λ[ind,:])
    end
    return p
end
