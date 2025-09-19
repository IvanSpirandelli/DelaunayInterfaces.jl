using GLMakie
using GLMakie.Makie.GeometryBasics
using GLMakie.Colors

fig = Figure()

interface_colormap = :viridis

gl_A = GridLayout(fig[1,1])
scene_A = LScene(gl_A[1, 1], show_axis = false)
ex_A = TopologicalInterfaces.EXAMPLES[:tetra_3_1]
TopologicalInterfaces._populate_interface_scene!(scene_A, ex_A.points, ex_A.labels, Vector{Float64}(); 
show_wireframe=true, show_mc_edges=true, show_mc_points=true, interface_colormap=interface_colormap)

gl_B = GridLayout(fig[1,2])
scene_B = LScene(gl_B[1, 1], show_axis = false)
ex_B = TopologicalInterfaces.EXAMPLES[:tetra_2_2]
TopologicalInterfaces._populate_interface_scene!(scene_B, ex_B.points, ex_B.labels, Vector{Float64}(); 
show_wireframe=true, show_mc_edges=true, show_mc_points=true, interface_colormap=interface_colormap)

gl_C = GridLayout(fig[2,1])
scene_C = LScene(gl_C[1, 1], show_axis = false)
ex_C = TopologicalInterfaces.EXAMPLES[:tetra_2_1_1]
TopologicalInterfaces._populate_interface_scene!(scene_C, ex_C.points, ex_C.labels, Vector{Float64}(); 
show_wireframe=true, show_mc_edges=true, show_mc_points=true, interface_colormap=interface_colormap)

gl_D = GridLayout(fig[2,2])
scene_D = LScene(gl_D[1, 1], show_axis = false)
ex_D = TopologicalInterfaces.EXAMPLES[:tetra_1_1_1_1]
TopologicalInterfaces._populate_interface_scene!(scene_D, ex_D.points, ex_D.labels, Vector{Float64}(); 
show_wireframe=true, show_mc_edges=true, show_mc_points=true, interface_colormap=interface_colormap)

for (label, layout) in zip(["A", "B", "C", "D"], [gl_A, gl_B, gl_C, gl_D])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 36,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
end

display(fig)