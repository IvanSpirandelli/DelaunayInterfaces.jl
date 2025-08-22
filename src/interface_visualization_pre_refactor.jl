using GLMakie
using GLMakie.Makie.GeometryBasics

# Your helper functions (get_barycentric_subdivision_and_filtration, get_multichromatic_tetrahedra, etc.)
# are assumed to be defined elsewhere in your code.

"""
    visualize_interface_sequence(
        points_sequence::Vector{Vector{Vector{Float64}}}, 
        labels_sequence::Vector{Vector{Int}}; 
        show_wireframe::Bool = false,
        show_points::Bool = true,
        show_edges::Bool = true,
        global_colorrange::Bool = true,
        colormap = :viridis
    )

Constructs a comprehensive, interactive figure for a sequence of interfaces with clean gradient edges.

This function pre-computes all visualization data for a smooth slider experience and includes several customizable options.

# Arguments
- `points_sequence`: A vector where each element is the list of points for one frame.
- `labels_sequence`: A vector where each element is the list of labels for one frame.

# Keyword Arguments
- `show_wireframe::Bool`: If `true`, overlays a white wireframe on the interface mesh.
- `show_points::Bool`: If `true`, displays the centered atom points.
- `show_edges::Bool`: If `true`, displays the multicolored edges with clean color gradients.
- `global_colorrange::Bool`: If `true`, the colormap is scaled globally across all frames. If `false`, it's scaled to each individual frame.
- `colormap`: The colormap to use for the interface mesh (e.g., `:viridis`, `:plasma`).
"""
function visualize_interface_sequence(
    points_sequence::Vector{Vector{Vector{Float64}}}, 
    labels_sequence::Vector{Vector{Int}}; 
    show_wireframe::Bool = false,
    show_points::Bool = true,
    show_edges::Bool = true,
    global_colorrange::Bool = true,
    colormap = :viridis
)
    n_interfaces = length(points_sequence)
    if n_interfaces != length(labels_sequence)
        error("The number of point sets must be equal to the number of label sets.")
    end

    # --- Pre-computation Step ---
    all_meshes = Vector{GeometryBasics.Mesh}(undef, n_interfaces)
    all_mesh_colors = Vector{Vector{Float64}}(undef, n_interfaces)
    all_individual_colorranges = Vector{Tuple{Float64, Float64}}(undef, n_interfaces)
    
    # Pre-computation storage for points and edges
    all_centered_points = show_points ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
    all_point_colors = show_points ? Vector{Vector{RGBf}}(undef, n_interfaces) : []
    all_edge_points = show_edges ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
    all_edge_colors = show_edges ? Vector{Vector{RGBf}}(undef, n_interfaces) : []

    println("Pre-computing visualization data for $n_interfaces interfaces...")

    # Define the exact color mapping for points and edges from your example
    point_edge_color_map = let
        cm = cgrad(:Accent_6, 6, categorical=true)
        Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
    end

    for i in 1:n_interfaces
        points = points_sequence[i]
        labels = labels_sequence[i]

        # 1. Center points
        center = Point3f(sum(points) / length(points))
        centered_points = [Point3f(p) - center for p in points]
        
        if show_points
            all_centered_points[i] = centered_points
            all_point_colors[i] = [point_edge_color_map[l] for l in labels]
        end

        # 2. Compute interface mesh and colors
        bcs, filtration = get_barycentric_subdivision_and_filtration(points, labels)
        bcs_centered = [Point3f(bc) - center for bc in bcs]
        faces = [TriangleFace(e[1]) for e in filtration if length(e[1]) == 3]
        
        all_meshes[i] = GeometryBasics.Mesh(bcs_centered, faces)
        mesh_colors = [e[2] for e in filtration if length(e[1]) == 1]
        all_mesh_colors[i] = mesh_colors
        
        min_c, max_c = isempty(mesh_colors) ? (0.0, 1.0) : (minimum(mesh_colors), maximum(mesh_colors))
        all_individual_colorranges[i] = (min_c, max_c)

        # 3. Compute edge points and their explicit gradient colors
        if show_edges
            mc_tets = get_multichromatic_tetrahedra(points, labels)
            get_edges(tet) = [(tet[1], tet[2]), (tet[1], tet[3]), (tet[1], tet[4]), (tet[2], tet[3]), (tet[2], tet[4]), (tet[3], tet[4])]
            
            frame_edge_points = Point3f[]
            frame_edge_colors = RGBf[] # Store explicit colors, not labels
            
            for tet in eachrow(mc_tets)
                for (p1_idx, p2_idx) in get_edges(tet)
                    push!(frame_edge_points, centered_points[p1_idx], centered_points[p2_idx])
                    
                    # Get the exact start and end color for the gradient
                    color1 = point_edge_color_map[labels[p1_idx]]
                    color2 = point_edge_color_map[labels[p2_idx]]
                    push!(frame_edge_colors, color1, color2)
                end
            end
            all_edge_points[i] = frame_edge_points
            all_edge_colors[i] = frame_edge_colors
        end
    end
    println("Computation complete.")

    # --- Interactive Visualization Setup ---
    f = Figure(fontsize = 12)
    sl_i = Slider(f[2, 1], range = 1:n_interfaces, startvalue = 1)
    slider_value = sl_i.value
    cgl = GridLayout(f[1,1])

    # Observables for data selection
    current_mesh = lift(idx -> all_meshes[idx], slider_value)
    current_mesh_colors = lift(idx -> all_mesh_colors[idx], slider_value)

    final_colorrange_observable = if global_colorrange
        all_clr_values = vcat(all_mesh_colors...)
        (minimum(all_clr_values), maximum(all_clr_values))
    else
        lift(idx -> all_individual_colorranges[idx], slider_value)
    end

    i_sc = LScene(cgl[1, 1], show_axis=false, scenekw = (lights = [AmbientLight(RGBf(1.0, 1.0, 1.0))],))

    # Plot the interface mesh
    mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange_observable, colormap=colormap)
    if show_wireframe
        wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
    end

    # Plot edges if requested
    if show_edges
        current_edge_points = lift(idx -> all_edge_points[idx], slider_value)
        current_edge_colors = lift(idx -> all_edge_colors[idx], slider_value)
        # Plot segments with their explicit color vectors for perfect gradients
        linesegments!(i_sc, current_edge_points, color=current_edge_colors, linewidth=2)
    end

    # Plot points if requested
    if show_points
        current_points = lift(idx -> all_centered_points[idx], slider_value)
        current_point_colors = lift(idx -> all_point_colors[idx], slider_value)
        # Plot points with their explicit colors
        scatter!(i_sc, current_points, color=current_point_colors, markersize=15, strokewidth=1, strokecolor=:black)
    end
    
    f
end

function get_interface_visualization(points::Vector{Vector{Float64}}, n_atoms_per_mol::Int; show_wireframe = false)
    labels = get_labels(length(points), n_atoms_per_mol) 
    get_interface_visualization(points, labels, show_wireframe = show_wireframe)
end

function get_interface_visualization(points::Vector{Vector{Float64}}, labels::Vector{Int}; show_wireframe = false)
    bcs, filtration = get_barycentric_subdivision_and_filtration(points, labels)

    max_v = maximum([v for (_, v) in filtration]) 
    min_v = minimum([v for (_, v) in filtration])
    
    f = Figure(fontsize = 12)

    i_sc = LScene(f[1:2, 1:2], show_axis=false, scenekw = (lights = [AmbientLight(RGBf(1.0, 1.0, 1.0))],))
    fs = [TriangleFace(e[1]) for e in filtration if length(e[1]) == 3]
    bcs = [Point3f(e) for e in bcs]
    msh = GeometryBasics.Mesh(bcs, fs)

    clr = [e[2] for e in filtration if length(e[1]) == 1]
    mesh!(i_sc, msh, color=clr, colorrange = (min_v, max_v), colormap = colormap)
    if show_wireframe
        wireframe!(i_sc, msh, color=:white, linewidth = 1)
    end
    f
end

function get_interface_and_multicolored_tetrahedron_visualization(points::Vector{Vector{Float64}}, n_atoms_per_mol::Int; show_mc_edges = false, show_wireframe = false, show_barycenters = false, interface_colormap = :viridis)
    labels = get_labels(length(points), n_atoms_per_mol) 
    get_interface_and_multicolored_tetrahedron_visualization(points, labels, show_mc_edges = show_mc_edges, show_wireframe = show_wireframe, show_barycenters = show_barycenters, interface_colormap = interface_colormap)
end

function get_interface_and_multicolored_tetrahedron_visualization(points::Vector{Vector{Float64}}, labels::Vector{Int}; show_mc_edges = false, show_wireframe = false, show_barycenters = false, interface_colormap = :viridis)
    bcs, filtration = get_barycentric_subdivision_and_filtration(points, labels)
    get_edges(tet) = [
        (tet[1], tet[2]),
        (tet[1], tet[3]),
        (tet[1], tet[4]),
        (tet[2], tet[3]),
        (tet[2], tet[4]),
        (tet[3], tet[4])
    ]

    f = Figure(fontsize = 12)
    ms = 15
    cm = :Dark2_8

    i_sc = LScene(f[1:2, 1:2], show_axis=false, scenekw = (lights = [AmbientLight(RGBf(1.0, 1.0, 1.0))],))
    fs = [TriangleFace(e[1]) for e in filtration if length(e[1]) == 3]

    bcs = [Point3f(e) for e in bcs]
    msh = GeometryBasics.Mesh(bcs, fs)
    clr = [e[2] for e in filtration if length(e[1]) == 1]

    max_v = maximum([c for c in clr])
    min_v = minimum([c for c in clr])

    mesh!(i_sc, msh, color=clr, colorrange = (min_v, max_v), colormap = interface_colormap)
    if show_wireframe
        wireframe!(i_sc, msh, color=:white, linewidth = 1)
    end

    if show_barycenters
        for (i, bc) in enumerate(bcs)
            scatter!(i_sc, [bc], color=:black, markersize = ms, overdraw = true)
            text!(i_sc, [bc]; text="$(i)", fontsize = 15, color=:red, overdraw = true)
        end
    end

    if show_mc_edges
        cm  = cgrad(:Accent_6, categorical=true)
        c1 = cm[1]
        c2 = cm[2]
        c3 = cm[3]
        c4 = cm[6]
        mc_tets = get_multichromatic_tetrahedra(points, labels)
        for tet in eachrow(mc_tets)
            for (i,j) in get_edges(tet)
                if i > j 
                    i, j = j, i
                end
                if Set([labels[i], labels[j]]) == Set([1,2])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c1, c2], linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([1,3])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c1, c3], linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([1,4])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c1, c4], linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([2,3])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c2, c3], linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([2,4])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c2, c4], linewidth=2)            
                elseif Set([labels[i], labels[j]]) == Set([3,4])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = 1:2, colormap=[c3, c4], linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([1])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = c1, linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([2])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = c2, linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([3])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = c3, linewidth=2)
                elseif Set([labels[i], labels[j]]) == Set([4])
                    lines!(i_sc, [Point3f(points[i]), Point3f(points[j])], color = c4, linewidth=2)
                end
            end
        end
    end
    f
end