points, radii, labels = TopologicalInterfaces.get_points_radii_labels("assets/pdbs/4bmg.pdb", ["A", "B"])
inflated_radii = (radii .+ 1.4) .* 1.1 # Inflate radii for visualization

molecules_colormap = cgrad(:Accent_6, 3, categorical=true)[[3,2]]
single_colormap = cgrad(:Accent_6, 3, categorical=true)[[2]]
interface_colormap = :viridis

fig = Figure()

gl_A = GridLayout(fig[1,1])
scene_A = LScene(gl_A[1, 1], show_axis = false)
TopologicalInterfaces._populate_molecules_scene!(scene_A, points, labels, radii; molecules_colormap=molecules_colormap)

gl_B = GridLayout(fig[2,1])
scene_B = LScene(gl_B[1, 1], show_axis = false)
TopologicalInterfaces._populate_interface_scene!(scene_B, points, labels, inflated_radii; interface_colormap=interface_colormap)

gl_C = GridLayout(fig[1:2,2:3])
scene_C = LScene(gl_C[1, 1], show_axis = false)
indices = findall(x -> x == 1, labels)
subset_points = points[indices]
subset_labels = labels[indices]
subset_radii = radii[indices] .* 0.75
TopologicalInterfaces._populate_molecules_scene!(scene_C, subset_points, subset_labels, subset_radii; molecules_colormap=single_colormap)
TopologicalInterfaces._populate_interface_scene!(scene_C, points, labels, inflated_radii; interface_colormap=interface_colormap)

for (label, layout) in zip(["A", "B", "C"], [gl_A, gl_B, gl_C])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 36,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
end

display(fig)