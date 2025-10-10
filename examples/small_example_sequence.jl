function display_example_sequence()
    points_sequence = [
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]],
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]],
        [[0.1, 0.1, 0.1], [2.0, 0.0, 0.2], [1.0, 2.5, 0.0], [1.0, 1.4, 1.0]],
        [[0.0, 0.0, 0.0], [2.0, 0.15, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 2.4]],
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2], [4.0, 2.0, 1.8], [2.5, 2.5, 2.5]],
    ]

    labels_sequence = [
        [1,1,1,2],
        [1,1,2,2],
        [1,1,2,3],
        [1,2,3,4],
        [1,1,1,2,2,2],
    ]

    #centered_points_sequence = [points .- sum(p -> p, points) / length(points) for points in points_sequence]

    interface_surfaces = [TopologicalInterfaces.InterfaceSurface(p,l) for (p,l) in zip(points_sequence, labels_sequence)]

    # Display the sequence visualization
    display(TopologicalInterfaces.visualize_interface_sequence(
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
