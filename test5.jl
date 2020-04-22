using Blink, MeshCat, GeometryTypes, CoordinateTransformations

vis = Visualizer()
open(vis, Blink.Window())
# open(vis)

setobject!(vis[:box1], 
    HyperRectangle(Vec(0., 0, 0), Vec(0.1, 0.2, 0.3)))

anim = Animation()

atframe(anim, 0) do
    settransform!(vis[:box1], Translation(0., 0, 0))
end
atframe(anim, 30) do
    settransform!(vis[:box1], Translation(0., 1, 0))
end

setanimation!(vis, anim)
