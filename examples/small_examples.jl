using GLMakie

SMALL_EXAMPLES = Dict{String, Any}(
    "tetra_3_1" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]], "color_labels" => [1,1,1,2]),
    "tetra_2_2" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]], "color_labels" => [1,1,2,2]),
    "tetra_2_1_1" => Dict("points" => [[0.1, 0.1, 0.1], [2.0, 0.0, 0.2], [1.0, 2.5, 0.0], [1.0, 1.4, 1.0]], "color_labels" => [1,1,2,3]),
    "tetra_1_1_1_1" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.15, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 2.4]], "color_labels" => [1,2,3,4]),
    "sheet_3_3" => Dict("points" => [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2], [4.0, 2.0, 1.8], [2.5, 2.5, 2.5]], "color_labels" => [1,1,1,2,2,2])
)

function display_tetrahedron(example_id = "tetra_1_1_1_1")
    points = SMALL_EXAMPLES[example_id]["points"]
    color_labels = SMALL_EXAMPLES[example_id]["color_labels"]

    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels)
    display(TopologicalInterfaces.get_interface_figure(interface_surface, points, color_labels; show_multicolored_points = true, show_multicolored_edges = true))
end

function display_all_tetrahedra()
    fig = Figure()

    interface_colormap = :viridis

    gl_A = GridLayout(fig[1, 1])
    scene_A = LScene(gl_A[1, 1], show_axis = false)
    points_A = SMALL_EXAMPLES["tetra_3_1"]["points"]
    color_labels_A = SMALL_EXAMPLES["tetra_3_1"]["color_labels"]
    interface_A = TopologicalInterfaces.InterfaceSurface(points_A, color_labels_A)
    TopologicalInterfaces.draw_interface_to_scene!(scene_A, interface_A, points_A, color_labels_A, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_B = GridLayout(fig[1, 2])
    scene_B = LScene(gl_B[1, 1], show_axis = false)
    points_B = SMALL_EXAMPLES["tetra_2_2"]["points"]
    color_labels_B = SMALL_EXAMPLES["tetra_2_2"]["color_labels"]
    interface_B = TopologicalInterfaces.InterfaceSurface(points_B, color_labels_B)
    TopologicalInterfaces.draw_interface_to_scene!(scene_B, interface_B, points_B, color_labels_B, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_C = GridLayout(fig[2, 1])
    scene_C = LScene(gl_C[1, 1], show_axis = false)    
    points_C = SMALL_EXAMPLES["tetra_2_1_1"]["points"]
    color_labels_C = SMALL_EXAMPLES["tetra_2_1_1"]["color_labels"]
    interface_C = TopologicalInterfaces.InterfaceSurface(points_C, color_labels_C)
    TopologicalInterfaces.draw_interface_to_scene!(scene_C, interface_C, points_C, color_labels_C, Vector{Float64}(); 
    show_wireframe=true, show_multicolored_edges=true, show_multicolored_points=true, interface_colormap=interface_colormap)

    gl_D = GridLayout(fig[2, 2])
    scene_D = LScene(gl_D[1, 1], show_axis = false)
    points_D = SMALL_EXAMPLES["tetra_1_1_1_1"]["points"]
    color_labels_D = SMALL_EXAMPLES["tetra_1_1_1_1"]["color_labels"]
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

function display_example_sequence(keys = ["tetra_3_1", "tetra_2_2", "tetra_2_1_1", "tetra_1_1_1_1", "sheet_3_3"])
    points_sequence = [SMALL_EXAMPLES[k]["points"] for k in keys]
    labels_sequence = [SMALL_EXAMPLES[k]["color_labels"] for k in keys]

    #centered_points_sequence = [points .- sum(p -> p, points) / length(points) for points in points_sequence]

    interface_surfaces = [TopologicalInterfaces.InterfaceSurface(p,l) for (p,l) in zip(points_sequence, labels_sequence)]

    display(TopologicalInterfaces.get_interface_sequence_figure(
        interface_surfaces,
        points_sequence,
        labels_sequence,
        fill(Vector{Float64}(), length(points_sequence));
        show_wireframe = true,
        show_multicolored_points = true, 
        show_multicolored_edges = true,
        global_colorrange = false
    ))
end
