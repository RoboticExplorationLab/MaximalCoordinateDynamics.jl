mutable struct Translational3{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    cid::Int64

    function Translational3(body1::AbstractBody{T}, body2::AbstractBody{T}, p1::AbstractVector{T}, p2::AbstractVector{T}) where T
        Nc = 3
        vertices = (p1, p2)
        cid = body2.id

        new{T,Nc}(vertices, cid), body1.id, body2.id
    end
end

function setForce!(joint::Translational3, body1::AbstractBody, body2::Body, F, No)
    return
end

@inline function minimalCoordinates(joint::Translational3, body1::AbstractBody{T}, body2::Body, No) where T
    SVector{0,T}()
end

@inline function g(joint::Translational3, body1::Body, body2::Body, Δt, No)
    vertices = joint.vertices
    getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - (getx3(body1, Δt) + vrotate(vertices[1], getq3(body1, Δt)))
end

@inline function ∂g∂posa(joint::Translational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = SMatrix{3,3,T,9}(-I)

        q = body1.q[No]
        R = -2 * VRᵀmat(q) * Rmat(Quaternion(joint.vertices[1])) * LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational3{T}, body1::AbstractBody, body2::Body, No) where T
    if body2.id == joint.cid
        X = SMatrix{3,3,T,9}(I)

        q = body2.q[No]
        R = 2 * VRᵀmat(q) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = SMatrix{3,3,T,9}(-Δt * I)

        q = body1.q[No]
        Ω = -2 * VRᵀmat(q) * Lmat(q) * Rᵀmat(ωbar(body1, Δt)) * Rmat(Quaternion(joint.vertices[1])) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3{T}, body1::AbstractBody, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = SMatrix{3,3,T,9}(Δt * I)

        q = body2.q[No]
        Ω = 2 * VRᵀmat(q) * Lmat(q) * Rᵀmat(ωbar(body2, Δt)) * Rmat(Quaternion(joint.vertices[2])) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function minimalCoordinates(joint::Translational3, body1::Origin{T}, body2::Body, No) where T
    SVector{0,T}()
end

@inline function g(joint::Translational3, body1::Origin, body2::Body, Δt, No)
    vertices = joint.vertices
    getx3(body2, Δt) + vrotate(vertices[2], getq3(body2, Δt)) - vertices[1]
end