using GLMakie
using GLMakie.Makie.GeometryBasics
using GLMakie.Colors

TETRAHEDRA = Dict{String, Any}(
    "tetra_3_1" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]], "color_labels" => [1,1,1,2]),
    "tetra_2_2" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]], "color_labels" => [1,1,2,2]),
    "tetra_2_1_1" => Dict("points" => [[0.1, 0.1, 0.1], [2.0, 0.0, 0.2], [1.0, 2.5, 0.0], [1.0, 1.4, 1.0]], "color_labels" => [1,1,2,3]),
    "tetra_1_1_1_1" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.15, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 2.4]], "color_labels" => [1,2,3,4]),
)

function display_tetrahedron(example_id = "tetra_1_1_1_1")
    points = TETRAHEDRA[example_id]["points"]
    color_labels = TETRAHEDRA[example_id]["color_labels"]

    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels)
    display(TopologicalInterfaces.get_interface_figure(interface_surface, points, color_labels; show_multicolored_points = true, show_multicolored_edges = true))
end

function display_all_tetrahedra_variants()
    fig = Figure()

    interface_colormap = :viridis

    gl_A = GridLayout(fig[1, 1])
    scene_A = LScene(gl_A[1, 1], show_axis = false)
    points_A = TETRAHEDRA["tetra_3_1"]["points"]
    color_labels_A = TETRAHEDRA["tetra_3_1"]["color_labels"]
    interface_A = TopologicalInterfaces.InterfaceSurface(points_A, color_labels_A)
    TopologicalInterfaces.draw_interface_to_scene!(scene_A, interface_A, points_A, color_labels_A, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_B = GridLayout(fig[1, 2])
    scene_B = LScene(gl_B[1, 1], show_axis = false)
    points_B = TETRAHEDRA["tetra_2_2"]["points"]
    color_labels_B = TETRAHEDRA["tetra_2_2"]["color_labels"]
    interface_B = TopologicalInterfaces.InterfaceSurface(points_B, color_labels_B)
    TopologicalInterfaces.draw_interface_to_scene!(scene_B, interface_B, points_B, color_labels_B, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_C = GridLayout(fig[2, 1])
    scene_C = LScene(gl_C[1, 1], show_axis = false)    
    points_C = TETRAHEDRA["tetra_2_1_1"]["points"]
    color_labels_C = TETRAHEDRA["tetra_2_1_1"]["color_labels"]
    interface_C = TopologicalInterfaces.InterfaceSurface(points_C, color_labels_C)
    TopologicalInterfaces.draw_interface_to_scene!(scene_C, interface_C, points_C, color_labels_C, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_D = GridLayout(fig[2, 2])
    scene_D = LScene(gl_D[1, 1], show_axis = false)
    points_D = TETRAHEDRA["tetra_1_1_1_1"]["points"]
    color_labels_D = TETRAHEDRA["tetra_1_1_1_1"]["color_labels"]
    interface_D = TopologicalInterfaces.InterfaceSurface(points_D, color_labels_D)
    TopologicalInterfaces.draw_interface_to_scene!(scene_D, interface_D, points_D, color_labels_D, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    for (label, layout) in zip(["A", "B", "C", "D"], [gl_A, gl_B, gl_C, gl_D])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 36,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    display(fig)
end