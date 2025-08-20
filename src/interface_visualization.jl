using GLMakie
using GLMakie.Makie.GeometryBasics

#================================================================================#
#                         INTERNAL HELPER FUNCTIONS                              #
#================================================================================#

"""
    _compute_and_prepare_mesh_data(points, labels; center_point=nothing)

Internal function to compute the barycentric subdivision and generate a renderable mesh.

It optionally centers the mesh vertices around a given `center_point`.

# Arguments
- `points`: The list of atom coordinates.
- `labels`: The list of atom labels.
- `center_point`: An optional `Point3f` to be subtracted from each mesh vertex.

# Returns
- `GeometryBasics.Mesh`: The renderable interface mesh.
- `Vector{Float64}`: A vector of scalar values for coloring the mesh vertices.
- `Vector{Point3f}`: The coordinates of the barycenter vertices (potentially centered).
"""
function _compute_and_prepare_mesh_data(points::Vector{Vector{Float64}}, labels::Vector{Int}; center_point=nothing)
    bcs, filtration = get_barycentric_subdivision_and_filtration(points, labels)
    
    # Extract faces (triangles) and vertex colors from the filtration data
    faces = [TriangleFace(e[1]) for e in filtration if length(e[1]) == 3]
    mesh_colors = [e[2] for e in filtration if length(e[1]) == 1]
    
    # Convert barycentric coordinates to Point3f
    bcs_points = [Point3f(bc) for bc in bcs]
    
    # Center the points if a center_point is provided
    if !isnothing(center_point)
        bcs_points .-= center_point
    end
    
    mesh = GeometryBasics.Mesh(bcs_points, faces)
    
    return mesh, mesh_colors, bcs_points
end

"""
    _prepare_mc_edge_visuals(points, labels, mc_tets, color_map; center_point=nothing)

Internal function to prepare vertex and color data for rendering multicolored tetrahedron edges.

This function creates two corresponding vectors: one for line segment endpoints and one for their
colors, allowing for clean gradients to be rendered efficiently.

# Arguments
- `points`: The list of original atom coordinates.
- `labels`: The list of atom labels.
- `mc_tets`: The multicolored tetrahedra data.
- `color_map`: A dictionary mapping labels to specific colors.
- `center_point`: An optional `Point3f` to be subtracted from each point.

# Returns
- `Vector{Point3f}`: A flat vector of points for `linesegments!`.
- `Vector{RGBf}`: A flat vector of colors corresponding to each point.
"""
function _prepare_mc_edge_visuals(points::Vector{Vector{Float64}}, labels::Vector{Int}, mc_tets, color_map; center_point=nothing)
    # Helper to get all 6 edges from a tetrahedron's vertex indices
    get_edges(tet) = [(tet[1], tet[2]), (tet[1], tet[3]), (tet[1], tet[4]), 
                      (tet[2], tet[3]), (tet[2], tet[4]), (tet[3], tet[4])]

    # Convert original points to Point3f and center if needed
    points3f = [Point3f(p) for p in points]
    if !isnothing(center_point)
        points3f .-= center_point
    end

    edge_points = Point3f[]
    edge_colors = RGBf[]

    for tet in eachrow(mc_tets)
        for (p1_idx, p2_idx) in get_edges(tet)
            # Add the start and end points of the edge
            push!(edge_points, points3f[p1_idx], points3f[p2_idx])
            
            # Add the start and end colors for the gradient
            color1 = color_map[labels[p1_idx]]
            color2 = color_map[labels[p2_idx]]
            push!(edge_colors, color1, color2)
        end
    end
    
    return edge_points, edge_colors
end

"""
    _plot_barycenters!(scene, bcs_points; markersize=15)

Internal plotting function to draw barycenters as scatter points with index labels.
"""
function _plot_barycenters!(scene, bcs_points; markersize=15)
    scatter!(scene, bcs_points, color=:black, markersize=markersize, overdraw=true)
    # Add text labels for each barycenter index
    for (i, bc) in enumerate(bcs_points)
        text!(scene, bc; text="$(i)", fontsize=15, color=:red, overdraw=true, align=(:center, :center))
    end
end


#================================================================================#
#                           PUBLIC API FUNCTIONS                                 #
#================================================================================#

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
    all_centered_points = show_points ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
    all_point_colors = show_points ? Vector{Vector{RGBf}}(undef, n_interfaces) : []
    all_edge_points = show_edges ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
    all_edge_colors = show_edges ? Vector{Vector{RGBf}}(undef, n_interfaces) : []

    println("Pre-computing visualization data for $n_interfaces interfaces...")

    point_edge_color_map = let
        cm = cgrad(:Accent_6, 6, categorical=true)
        Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
    end

    for i in 1:n_interfaces
        points = points_sequence[i]
        labels = labels_sequence[i]

        # 1. Center points for consistent visualization
        center = Point3f(sum(p -> Point3f(p), points) / length(points))
        
        if show_points
            all_centered_points[i] = [Point3f(p) - center for p in points]
            all_point_colors[i] = [point_edge_color_map[l] for l in labels]
        end

        # 2. Compute interface mesh using the helper function
        mesh, mesh_colors, _ = _compute_and_prepare_mesh_data(points, labels, center_point=center)
        all_meshes[i] = mesh
        all_mesh_colors[i] = mesh_colors
        
        min_c, max_c = isempty(mesh_colors) ? (0.0, 1.0) : (minimum(mesh_colors), maximum(mesh_colors))
        all_individual_colorranges[i] = (min_c, max_c)

        # 3. Compute edge data using the helper function
        if show_edges
            mc_tets = get_multichromatic_tetrahedra(points, labels)
            edge_points, edge_colors = _prepare_mc_edge_visuals(points, labels, mc_tets, point_edge_color_map, center_point=center)
            all_edge_points[i] = edge_points
            all_edge_colors[i] = edge_colors
        end
    end
    println("Computation complete.")

    # --- Interactive Visualization Setup ---
    f = Figure(fontsize=12)
    sl_i = Slider(f[2, 1], range=1:n_interfaces, startvalue=1)
    slider_value = sl_i.value
    cgl = GridLayout(f[1,1])
    i_sc = LScene(cgl[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

    # Observables for interactive data
    current_mesh = lift(idx -> all_meshes[idx], slider_value)
    current_mesh_colors = lift(idx -> all_mesh_colors[idx], slider_value)
    final_colorrange = if global_colorrange
        (minimum(vcat(all_mesh_colors...)), maximum(vcat(all_mesh_colors...)))
    else
        lift(idx -> all_individual_colorranges[idx], slider_value)
    end

    # Plot interface mesh
    mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange, colormap=colormap)
    if show_wireframe
        wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
    end

    # Plot edges if requested
    if show_edges
        current_edge_points = lift(idx -> all_edge_points[idx], slider_value)
        current_edge_colors = lift(idx -> all_edge_colors[idx], slider_value)
        linesegments!(i_sc, current_edge_points, color=current_edge_colors, linewidth=2)
    end

    # Plot points if requested
    if show_points
        current_points = lift(idx -> all_centered_points[idx], slider_value)
        current_point_colors = lift(idx -> all_point_colors[idx], slider_value)
        scatter!(i_sc, current_points, color=current_point_colors, markersize=15, strokewidth=1, strokecolor=:black)
    end
    
    return f
end

"""
Visualizes a single interface, deriving labels from the number of atoms per molecule.
"""
function get_interface_visualization(points::Vector{Vector{Float64}}, n_atoms_per_mol::Int; show_wireframe = false, colormap = :viridis)
    labels = get_labels(length(points), n_atoms_per_mol) 
    return get_interface_visualization(points, labels, show_wireframe=show_wireframe, colormap=colormap)
end

"""
Visualizes a single, static interface mesh.
"""
function get_interface_visualization(points::Vector{Vector{Float64}}, labels::Vector{Int}; show_wireframe = false, colormap = :viridis)
    # 1. Compute mesh data (no centering for this static view)
    msh, clr, _ = _compute_and_prepare_mesh_data(points, labels)

    # 2. Determine color range for the mesh
    min_v, max_v = isempty(clr) ? (0.0, 1.0) : (minimum(clr), maximum(clr))
    
    # 3. Setup figure and scene
    f = Figure(fontsize=12)
    i_sc = LScene(f[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))
    
    # 4. Plot mesh and optional wireframe
    mesh!(i_sc, msh, color=clr, colorrange=(min_v, max_v), colormap=colormap)
    if show_wireframe
        wireframe!(i_sc, msh, color=:white, linewidth=1)
    end
    
    return f
end

"""
Visualizes an interface with optional multicolored tetrahedra, deriving labels automatically.
"""
function get_interface_and_multicolored_tetrahedron_visualization(points::Vector{Vector{Float64}}, n_atoms_per_mol::Int; show_mc_edges = false, show_wireframe = false, show_barycenters = false, interface_colormap = :viridis)
    labels = get_labels(length(points), n_atoms_per_mol) 
    return get_interface_and_multicolored_tetrahedron_visualization(points, labels, show_mc_edges=show_mc_edges, show_wireframe=show_wireframe, show_barycenters=show_barycenters, interface_colormap=interface_colormap)
end

"""
Visualizes a single, static interface with options to show multicolored tetrahedra and barycenters.
"""
function get_interface_and_multicolored_tetrahedron_visualization(points::Vector{Vector{Float64}}, labels::Vector{Int}; show_mc_edges = false, show_wireframe = false, show_barycenters = false, interface_colormap = :viridis)
    # 1. Setup figure and scene
    f = Figure(fontsize=12)
    i_sc = LScene(f[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

    # 2. Compute and plot the interface mesh
    msh, clr, bcs_points = _compute_and_prepare_mesh_data(points, labels)
    min_v, max_v = isempty(clr) ? (0.0, 1.0) : (minimum(clr), maximum(clr))
    
    mesh!(i_sc, msh, color=clr, colorrange=(min_v, max_v), colormap=interface_colormap)
    if show_wireframe
        wireframe!(i_sc, msh, color=:white, linewidth=1)
    end

    # 3. Plot barycenters if requested
    if show_barycenters
        _plot_barycenters!(i_sc, bcs_points, markersize=15)
    end

    # 4. Compute and plot multicolored edges if requested
    if show_mc_edges
        color_map = let
            cm = cgrad(:Accent_6, 6, categorical=true)
            Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
        end
        
        mc_tets = get_multichromatic_tetrahedra(points, labels)
        edge_points, edge_colors = _prepare_mc_edge_visuals(points, labels, mc_tets, color_map)
        
        linesegments!(i_sc, edge_points, color=edge_colors, linewidth=2)
    end
    
    return f
end