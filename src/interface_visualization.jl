using GLMakie
using GLMakie.Makie.GeometryBasics
using GLMakie.Colors

function get_interface_and_point_cloud_figure(
    interface_surface::InterfaceSurface,
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64} = Vector{Float64}(); 
    show_axis = false, 
    show_wireframe = false, 
    interface_colormap = :viridis,
    point_cloud_colormap = :Dark2_4
    )
    f = Figure()


    point_cloud_gl = GridLayout(f[1,1])
    point_cloud_lscene = LScene(point_cloud_gl[1,1], show_axis = show_axis)
    draw_point_cloud_to_scene!(point_cloud_lscene, points, color_labels, radii; point_cloud_colormap)

    interface_gl = GridLayout(f[1,2])
    interface_lscene = LScene(interface_gl[1,1], show_axis = show_axis)
    draw_interface_to_scene!(interface_lscene, interface_surface; show_wireframe, interface_colormap)

    for (label, layout) in zip(["A", "B"], [interface_gl, point_cloud_gl])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 36,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end
    f
end

function get_interface_figure(
    interface_surface::InterfaceSurface; 
    show_axis = false, 
    show_wireframe = false, 
    show_barycenters = false, 
    interface_colormap = :viridis,
    )
    f = Figure()
    lscene = LScene(f[1,1], show_axis = show_axis)
    draw_interface_to_scene!(lscene, interface_surface; show_wireframe, show_barycenters, interface_colormap)
    _draw_multicolored_points_to_scene!
    f
end

function get_interface_figure(
    interface_surface::InterfaceSurface,
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64} = Vector{Float64}(); 
    show_axis = false, 
    show_wireframe = false, 
    show_barycenters = false, 
    interface_colormap = :viridis,
    show_multicolored_points = false,
    show_multicolored_edges = false,
    )
    f = Figure()
    lscene = LScene(f[1,1], show_axis = show_axis)
    draw_interface_to_scene!(
        lscene, 
        interface_surface,
        points,
        color_labels,
        radii; 
        show_wireframe, 
        show_barycenters, 
        interface_colormap,
        show_multicolored_points,
        show_multicolored_edges,
    )
    f
end

function draw_interface_to_scene!(
    lscene::LScene, 
    interface_surface::InterfaceSurface;
    show_wireframe = false, 
    show_barycenters = false,
    interface_colormap = :viridis
)
    mesh, mesh_colors = _generate_colored_mesh(interface_surface)

    min_v, max_v = (minimum(mesh_colors), maximum(mesh_colors))
    
    mesh!(lscene, mesh, color=mesh_colors, colorrange=(min_v, max_v), colormap=interface_colormap, shading = false)
    if show_wireframe
        wireframe!(lscene, mesh, color=:white, linewidth=1)
    end

    if show_barycenters
        _draw_barycenters_to_scene!(lscene, [Point3f(v) for v in interface_surface.vertices])
    end
end

function draw_interface_to_scene!(
    lscene::LScene, 
    interface_surface::InterfaceSurface,
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64} = Vector{Float64}(); 
    show_wireframe = false, 
    show_barycenters = false, 
    interface_colormap = :viridis,
    show_multicolored_points = false,
    show_multicolored_edges = false,
)
    mesh, mesh_colors = _generate_colored_mesh(interface_surface)

    min_v, max_v = (minimum(mesh_colors), maximum(mesh_colors))
    
    mesh!(lscene, mesh, color=mesh_colors, colorrange=(min_v, max_v), colormap=interface_colormap, shading = false)
    if show_wireframe
        wireframe!(lscene, mesh, color=:white, linewidth=1)
    end

    if show_barycenters
        _draw_barycenters_to_scene!(lscene, [Point3f(v) for v in interface_surface.vertices])
    end

    if show_multicolored_points
        _draw_multicolored_points_to_scene!(lscene, points, color_labels)
    end

    if show_multicolored_edges 
        _draw_multicolored_edges_to_scene!(lscene, points, color_labels, radii, interface_surface.weighted, interface_surface.alpha)
    end
end

function get_point_cloud_figure(
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64} = Vector{Float64}(); 
    show_axis = false, 
    point_cloud_colormap = :Dark2_4
    )
    f = Figure()
    lscene = LScene(f[1,1], show_axis = show_axis)
    draw_point_cloud_to_scene!(lscene, points, color_labels, radii; point_cloud_colormap = point_cloud_colormap)
    _draw_multicolored_points_to_scene!
    f
end

function draw_point_cloud_to_scene!(
    lscene::LScene,
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64} = Vector{Float64}();
    point_cloud_colormap = point_cloud_colormap
)
    meshscatter!(lscene, 
        [Point3f(p) for p in points], 
        markersize = radii,
        color=color_labels,
        colormap = point_cloud_colormap,
    ) 
end

function _draw_multicolored_points_to_scene!(
    lscene::LScene, 
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int};
    point_cloud_colormap = :Dark2_4
)
    scatter!(lscene, [Point3f(p) for p in points], color=cgrad(point_cloud_colormap, 4, categorical=true)[color_labels], strokecolor=:black, strokewidth=1)
end

function _draw_multicolored_edges_to_scene!(
    lscene::LScene, 
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64}, 
    weighted::Bool, 
    alpha::Bool;
    point_cloud_colormap = :Dark2_4
    )
    edge_points, edge_colors = _get_multicolored_edge_data(points, color_labels, radii, weighted, alpha; point_cloud_colormap=point_cloud_colormap)
    linesegments!(lscene, edge_points, color=edge_colors, linewidth=2)
end

function _get_multicolored_edge_data(
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64}, 
    weighted::Bool, 
    alpha::Bool;
    point_cloud_colormap = :Dark2_4
)
    mc_tets = get_multicolored_tetrahedra(points, color_labels, radii, weighted, alpha)
    # Helper to get all 6 edges from a tetrahedron's vertex indices
    get_edges(tet) = [(tet[1], tet[2]), (tet[1], tet[3]), (tet[1], tet[4]), 
                      (tet[2], tet[3]), (tet[2], tet[4]), (tet[3], tet[4])]

    # Convert original points to Point3f and center if needed
    points3f = [Point3f(p) for p in points]

    edge_points = Point3f[]
    edge_colors = RGBA[]

    color_map = let
        cm = cgrad(point_cloud_colormap, 4, categorical=true)
        Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[4])
    end

    for tet in mc_tets
        for (p1_idx, p2_idx) in get_edges(tet)
            # Add the start and end points of the edge
            push!(edge_points, points3f[p1_idx], points3f[p2_idx])
            
            # Add the start and end colors for the gradient
            color1 = color_map[color_labels[p1_idx]]
            color2 = color_map[color_labels[p2_idx]]
            push!(edge_colors, color1, color2)
        end
    end
    edge_points, edge_colors
end

function _draw_barycenters_to_scene!(scene, bcs_points; markersize=15)
    scatter!(scene, bcs_points, color=:black, markersize=markersize, overdraw=true)
    # Add text labels for each barycenter index
    for (i, bc) in enumerate(bcs_points)
        text!(scene, bc; text="$(i)", fontsize=15, color=:red, overdraw=true, align=(:center, :center))
    end
end

function _generate_colored_mesh(interface_surface::InterfaceSurface, max_value = Inf)
    # Extract faces (triangles) and vertex colors from the filtration data
    faces = [TriangleFace(e[1]) for e in interface_surface.filtration if length(e[1]) == 3 && e[2] <= max_value]
    mesh = GeometryBasics.Mesh([Point3f(v) for v in interface_surface.vertices], faces)

    mesh_colors = [e[2] for e in interface_surface.filtration if length(e[1]) == 1]
    return mesh, mesh_colors
end

function get_interface_filtration_figure(
    interface_surface::InterfaceSurface;
    show_wireframe::Bool = false,
    interface_surface_colormap = :viridis
)
    levels = sort!(unique([e[2] for e in interface_surface.filtration]))
    meshes_and_colors = [_generate_colored_mesh(interface_surface, lvl) for lvl in levels]
    all_meshes = [e[1] for e in meshes_and_colors]
    all_mesh_colors = [e[2] for e in meshes_and_colors]

    f = Figure(fontsize=12)
    sl_i = Slider(f[2, 1], range=1:length(levels), startvalue=length(levels))
    slider_value = sl_i.value
    cgl = GridLayout(f[1,1])
    i_sc = LScene(cgl[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

    # Observables for interactive data
    current_mesh = lift(idx -> all_meshes[idx], slider_value)
    current_mesh_colors = lift(idx -> all_mesh_colors[idx], slider_value)
    final_colorrange = (minimum(all_mesh_colors[end]), maximum(all_mesh_colors[end]))
    
    # Plot interface mesh
    mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange, colormap=interface_surface_colormap)
    if show_wireframe
        wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
    end
    f
end


function get_interface_sequence_figure(
    interface_surfaces::Vector{InterfaceSurface},
    points_sequence::Vector{Vector{Vector{Float64}}}, 
    color_labels_sequence::Vector{Vector{Int}},
    radii_sequence::Vector{Vector{Float64}};
    show_wireframe::Bool = false,
    show_multicolored_points::Bool = false,
    show_multicolored_edges::Bool = false,
    global_colorrange::Bool = false,
    interface_surface_colormap = :viridis,
    point_cloud_colormap = :Dark2_4
)
    n_interfaces = length(interface_surfaces)
    colored_meshes = [_generate_colored_mesh(is) for is in interface_surfaces]
    all_meshes = [m for (m,_) in colored_meshes]
    all_mesh_colors = [mc for (_,mc) in colored_meshes]

    all_individual_colorranges = [(minimum(mesh_colors), maximum(mesh_colors)) for mesh_colors in all_mesh_colors]

    all_edge_points = show_multicolored_edges ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
    all_edge_colors = show_multicolored_edges ? Vector{Vector{RGBf}}(undef, n_interfaces) : []
        
    if show_multicolored_edges
        for (i, ((points, labels), (radii, interface))) in enumerate(zip(zip(points_sequence, color_labels_sequence), zip(radii_sequence, interface_surfaces)))
            edge_points, edge_colors = _get_multicolored_edge_data(points, labels, radii, interface.weighted, interface.alpha)
            all_edge_points[i] = edge_points
            all_edge_colors[i] = edge_colors
        end
    end

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
    mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange, colormap=interface_surface_colormap)
    if show_wireframe
        wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
    end

    # Plot edges if requested
    if show_multicolored_edges
        current_edge_points = lift(idx -> all_edge_points[idx], slider_value)
        current_edge_colors = lift(idx -> all_edge_colors[idx], slider_value)
        linesegments!(i_sc, current_edge_points, color=current_edge_colors, linewidth=2)
    end

    # Plot points if requested
    if show_multicolored_points
        current_points = lift(idx -> [Point3f(e) for e in points_sequence[idx]], slider_value)
        current_point_colors = lift(idx -> cgrad(point_cloud_colormap, 4, categorical=true)[color_labels_sequence[idx]], slider_value)
        scatter!(i_sc, current_points, color=current_point_colors, markersize=15, strokewidth=1, strokecolor=:black)
    end
    
    return f
end
