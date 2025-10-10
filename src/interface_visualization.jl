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
    point_cloud_colormap = :Accent_6
    )
    f = Figure()
    interface_gl = GridLayout(f[1,1])
    interface_lscene = LScene(interface_gl[1,1], show_axis = show_axis)
    draw_interface_to_scene!(interface_lscene, interface_surface; show_wireframe, interface_colormap)

    point_cloud_gl = GridLayout(f[1,2])
    point_cloud_lscene = LScene(point_cloud_gl[1,1], show_axis = show_axis)
    draw_point_cloud_to_scene!(point_cloud_lscene, points, color_labels, radii; point_cloud_colormap)
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
    point_cloud_colormap = :Accent_6
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
        marker = Sphere(Point3f(0), 1),
        markersize = radii,
        color = color_labels,
        colormap = point_cloud_colormap
    ) 
end

function _draw_multicolored_points_to_scene!(
    lscene::LScene, 
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}
)
    color_map = let
        cm = cgrad(:Accent_6, 6, categorical=true)
        Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
    end
    
    scatter!(lscene, [Point3f(p) for p in points], color=[color_map[l] for l in color_labels], strokecolor=:black, strokewidth=1)
end

function _draw_multicolored_edges_to_scene!(
    lscene::LScene, 
    points::Vector{Vector{Float64}}, 
    color_labels::Vector{Int}, 
    radii::Vector{Float64}, 
    weighted::Bool, 
    alpha::Bool
    )
    color_map = let
        cm = cgrad(:Accent_6, 6, categorical=true)
        Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
    end

    mc_tets = get_multichromatic_tetrahedra(points, color_labels, radii, weighted, alpha)
    # Helper to get all 6 edges from a tetrahedron's vertex indices
    get_edges(tet) = [(tet[1], tet[2]), (tet[1], tet[3]), (tet[1], tet[4]), 
                      (tet[2], tet[3]), (tet[2], tet[4]), (tet[3], tet[4])]

    # Convert original points to Point3f and center if needed
    points3f = [Point3f(p) for p in points]

    edge_points = Point3f[]
    edge_colors = RGBf[]

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
    
    linesegments!(lscene, edge_points, color=edge_colors, linewidth=2)
end

function _draw_barycenters_to_scene!(scene, bcs_points; markersize=15)
    scatter!(scene, bcs_points, color=:black, markersize=markersize, overdraw=true)
    # Add text labels for each barycenter index
    for (i, bc) in enumerate(bcs_points)
        text!(scene, bc; text="$(i)", fontsize=15, color=:red, overdraw=true, align=(:center, :center))
    end
end

function _generate_colored_mesh(interface_surface::InterfaceSurface)
        # Extract faces (triangles) and vertex colors from the filtration data
        faces = [TriangleFace(e[1]) for e in interface_surface.filtration if length(e[1]) == 3]
        mesh = GeometryBasics.Mesh([Point3f(v) for v in interface_surface.vertices], faces)

        mesh_colors = [e[2] for e in interface_surface.filtration if length(e[1]) == 1]
        return mesh, mesh_colors
end


# function visualize_interface_sequence(
#     interface_surfaces::Vector{InterfaceSurface},
#     points_sequence::Vector{Vector{Vector{Float64}}}, 
#     labels_sequence::Vector{Vector{Int}};
#     show_wireframe::Bool = false,
#     show_points::Bool = false,
#     show_edges::Bool = false,
#     global_colorrange::Bool = false,
#     interface_surface_colormap = :viridis
# )
#     n_interfaces = length(interface_surfaces)
#     colored_meshes = [_generate_colored_mesh(is) for is in interface_surfaces]
#     all_meshes = [m for (m,_) in colored_meshes]
#     all_mesh_colors = [mc for (_,mc) in colored_meshes]

#     all_individual_colorranges = [(minimum(mesh_colors), maximum(mesh_colors)) for mesh_colors in all_mesh_colors]

#     all_point_colors = show_points ? Vector{Vector{RGBf}}(undef, n_interfaces) : []
#     all_edge_points = show_edges ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
#     all_edge_colors = show_edges ? Vector{Vector{RGBf}}(undef, n_interfaces) : []

#     point_edge_color_map = let
#         cm = cgrad(:Accent_6, 6, categorical=true)
#         Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
#     end

#     if show_points
#         all_point_colors = [[point_edge_color_map[l] for l in labels] for labels in labels_sequence]
#     end
        
#     if show_edges
#         for (i, ((points, labels), (interfaces))) in enumerate(zip(zip(points_sequence, labels_sequence), interface_surfaces))
#             mc_tets = get_multichromatic_tetrahedra(points, labels, radii, )
#             edge_points, edge_colors = _prepare_mc_edge_visuals(points, labels, mc_tets, point_edge_color_map, center_point=center)
#             all_edge_points[i] = edge_points
#             all_edge_colors[i] = edge_colors
#         end
#     end

#     # --- Interactive Visualization Setup ---
#     f = Figure(fontsize=12)
#     sl_i = Slider(f[2, 1], range=1:n_interfaces, startvalue=1)
#     slider_value = sl_i.value
#     cgl = GridLayout(f[1,1])
#     i_sc = LScene(cgl[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

#     # Observables for interactive data
#     current_mesh = lift(idx -> all_meshes[idx], slider_value)
#     current_mesh_colors = lift(idx -> all_mesh_colors[idx], slider_value)
#     final_colorrange = if global_colorrange
#         (minimum(vcat(all_mesh_colors...)), maximum(vcat(all_mesh_colors...)))
#     else
#         lift(idx -> all_individual_colorranges[idx], slider_value)
#     end

#     # Plot interface mesh
#     mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange, colormap=interface_surface_colormap)
#     if show_wireframe
#         wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
#     end

#     # Plot edges if requested
#     if show_edges
#         current_edge_points = lift(idx -> all_edge_points[idx], slider_value)
#         current_edge_colors = lift(idx -> all_edge_colors[idx], slider_value)
#         linesegments!(i_sc, current_edge_points, color=current_edge_colors, linewidth=2)
#     end

#     # Plot points if requested
#     if show_points
#         current_points = lift(idx -> all_centered_points[idx], slider_value)
#         current_point_colors = lift(idx -> all_point_colors[idx], slider_value)
#         scatter!(i_sc, current_points, color=current_point_colors, markersize=15, strokewidth=1, strokecolor=:black)
#     end
    
#     return f
# end

# #================================================================================#
# #                           PUBLIC API FUNCTIONS                                 #
# #================================================================================#

# function visualize_interface_sequence(
#     points_sequence::Vector{Vector{Vector{Float64}}}, 
#     labels_sequence::Vector{Vector{Int}},
#     radii_sequence::Vector{Vector{Float64}} = Vector{Vector{Float64}}(); 
#     show_wireframe::Bool = false,
#     show_points::Bool = true,
#     show_edges::Bool = true,
#     global_colorrange::Bool = true,
#     colormap = :viridis
# )
#     n_interfaces = length(points_sequence)
#     if n_interfaces != length(labels_sequence)
#         error("The number of point sets must be equal to the number of label sets.")
#     end

#     # --- Pre-computation Step ---
#     all_meshes = Vector{GeometryBasics.Mesh}(undef, n_interfaces)
#     all_mesh_colors = Vector{Vector{Float64}}(undef, n_interfaces)
#     all_individual_colorranges = Vector{Tuple{Float64, Float64}}(undef, n_interfaces)
#     all_centered_points = show_points ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
#     all_point_colors = show_points ? Vector{Vector{RGBf}}(undef, n_interfaces) : []
#     all_edge_points = show_edges ? Vector{Vector{Point3f}}(undef, n_interfaces) : []
#     all_edge_colors = show_edges ? Vector{Vector{RGBf}}(undef, n_interfaces) : []

#     println("Pre-computing visualization data for $n_interfaces interfaces...")

#     point_edge_color_map = let
#         cm = cgrad(:Accent_6, 6, categorical=true)
#         Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
#     end

#     for i in 1:n_interfaces
#         points = points_sequence[i]
#         labels = labels_sequence[i]
#         radii = length(radii_sequence) > 0 ? radii_sequence[i] : Float64[]

#         # 1. Center points for consistent visualization
#         center = Point3f(sum(p -> Point3f(p), points) / length(points))
        
#         if show_points
#             all_centered_points[i] = [Point3f(p) - center for p in points]
#             all_point_colors[i] = [point_edge_color_map[l] for l in labels]
#         end

#         # 2. Compute interface mesh using the helper function
#         mesh, mesh_colors, _ = _compute_and_prepare_mesh_data(points, labels, radii, center_point=center)
#         all_meshes[i] = mesh
#         all_mesh_colors[i] = mesh_colors
        
#         min_c, max_c = isempty(mesh_colors) ? (0.0, 1.0) : (minimum(mesh_colors), maximum(mesh_colors))
#         all_individual_colorranges[i] = (min_c, max_c)

#         # 3. Compute edge data using the helper function
#         if show_edges
#             mc_tets = get_multichromatic_tetrahedra(points, labels)
#             edge_points, edge_colors = _prepare_mc_edge_visuals(points, labels, mc_tets, point_edge_color_map, center_point=center)
#             all_edge_points[i] = edge_points
#             all_edge_colors[i] = edge_colors
#         end
#     end
#     println("Computation complete.")

#     # --- Interactive Visualization Setup ---
#     f = Figure(fontsize=12)
#     sl_i = Slider(f[2, 1], range=1:n_interfaces, startvalue=1)
#     slider_value = sl_i.value
#     cgl = GridLayout(f[1,1])
#     i_sc = LScene(cgl[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

#     # Observables for interactive data
#     current_mesh = lift(idx -> all_meshes[idx], slider_value)
#     current_mesh_colors = lift(idx -> all_mesh_colors[idx], slider_value)
#     final_colorrange = if global_colorrange
#         (minimum(vcat(all_mesh_colors...)), maximum(vcat(all_mesh_colors...)))
#     else
#         lift(idx -> all_individual_colorranges[idx], slider_value)
#     end

#     # Plot interface mesh
#     mesh!(i_sc, current_mesh, color=current_mesh_colors, colorrange=final_colorrange, colormap=colormap)
#     if show_wireframe
#         wireframe!(i_sc, current_mesh, color=:white, linewidth=1)
#     end

#     # Plot edges if requested
#     if show_edges
#         current_edge_points = lift(idx -> all_edge_points[idx], slider_value)
#         current_edge_colors = lift(idx -> all_edge_colors[idx], slider_value)
#         linesegments!(i_sc, current_edge_points, color=current_edge_colors, linewidth=2)
#     end

#     # Plot points if requested
#     if show_points
#         current_points = lift(idx -> all_centered_points[idx], slider_value)
#         current_point_colors = lift(idx -> all_point_colors[idx], slider_value)
#         scatter!(i_sc, current_points, color=current_point_colors, markersize=15, strokewidth=1, strokecolor=:black)
#     end
    
#     return f
# end


# function visualize_interface(points::Vector{Vector{Float64}}, labels::Vector{Int}, radii::Vector{Float64} = Vector{Float64}(); 
#     show_wireframe = false, interface_colormap = :viridis, show_barycenters = true,
#     show_mc_edges = true)
#     f = Figure(fontsize=12)
#     i_sc = LScene(f[1, 1], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))
#     _populate_interface_scene!(
#         i_sc, 
#         points, 
#         labels, 
#         radii; 
#         show_wireframe=show_wireframe, 
#         interface_colormap=interface_colormap,
#         show_barycenters=show_barycenters,
#         show_mc_edges=show_mc_edges
#     )
#     f
# end

# function visualize_point_cloud(points::Vector{Vector{Float64}}, labels::Vector{Int}, radii::Vector{Float64} = Vector{Float64}(); point_cloud_colormap = :Accent_6)
#     f = Figure(fontsize=12)
#     i_sc = LScene(f[1, 1], show_axis=false,)
#     _populate_point_cloud_scene!(
#         i_sc, 
#         points, 
#         labels, 
#         radii; 
#         point_cloud_colormap=point_cloud_colormap
#     )
#     f
# end

# function visualize_overlay(
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int}, 
#     atom_radii::Vector{Float64},
#     sasa_radii::Vector{Float64};
#     point_cloud_colormap = :Blues_3,
#     interface_colormap = :viridis,
#     show_wireframe = false,
#     point_cloud_transparency = 0.025
# )
#     f = Figure(fontsize=12)

#     lscene = LScene(f[1, 1], show_axis=false)#, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

#     _populate_interface_scene!(
#         lscene, 
#         points, 
#         labels, 
#         sasa_radii; 
#         show_wireframe=show_wireframe, 
#         interface_colormap=interface_colormap
#     )

#     if !isempty(points) && (length(points) == length(atom_radii) == length(labels))
#         meshscatter!(lscene, 
#             [Point3f(p) for p in points], 
#             marker = Sphere(Point3f(0), 1),
#             markersize = atom_radii,
#             color = labels,
#             colormap = (point_cloud_colormap, point_cloud_transparency),
#             transparency = point_cloud_transparency < 1.0
#         )
#     end

#     return f
# end

# function visualize_interactive_toggle(
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int}, 
#     atom_radii::Vector{Float64},
#     sasa_radii::Vector{Float64};
#     point_cloud_colormap = :grays,
#     interface_colormap = :magma,
#     show_wireframe = false,
#     point_cloud_transparency = 1.0
# )
#     # --- 1. Setup Figure and Layout ---
#     f = Figure(fontsize=16)

#     lscene = LScene(f[1:4, 1:4], show_axis=false, scenekw = (camera = cam3d!,),)

#     # Create a grid layout for the UI controls on the right side
#     controls = f[1:4, 4] = GridLayout(tellwidth = false)
#     Label(controls[1, 1:2], "point_cloud", fontsize=18, font=:bold)

#     # --- 2. Plot the Static Interface ---
#     _populate_interface_scene!(
#         lscene, 
#         points, 
#         labels, 
#         sasa_radii; 
#         show_wireframe=show_wireframe, 
#         interface_colormap=interface_colormap
#     )

#     # --- 3. Group Data and Plot point_cloud Individually ---
#     unique_labels = sort(unique(labels))
#     point_cloud_plots = Dict{Int, Any}() # Dictionary to store plot objects for each label

#     categorical_colors = cgrad(point_cloud_colormap, length(unique_labels), categorical=true)

#     for (i, label) in enumerate(unique_labels)
#         # Find all indices corresponding to the current label
#         indices = findall(x -> x == label, labels)
        
#         # Skip if no points found for this label
#         if isempty(indices)
#             continue
#         end

#         # Extract the data for the current group
#         group_points = [Point3f(p) for p in points[indices]]
#         group_radii = atom_radii[indices]
        
#         # Plot this group and store the plot object in our dictionary
#         # We use the pre-calculated color for this group
#         point_cloud_plots[label] = if point_cloud_transparency >= 1.0
#             # OPAQUE rendering path
#             meshscatter!(lscene, group_points,
#                 marker=Sphere(Point3f(0), 1),
#                 markersize=group_radii,
#                 color=categorical_colors[i],
#                 transparency=false
#             )
#         else
#             # TRANSPARENT rendering path
#             meshscatter!(lscene, group_points,
#                 marker=Sphere(Point3f(0), 1),
#                 markersize=group_radii,
#                 color=(categorical_colors[i], point_cloud_transparency),
#                 transparency=true
#             )
#         end
#     end

#     # --- 4. Create UI Toggles and Link them to Plots ---
#     for (i, label) in enumerate(unique_labels)
#         # Create a toggle and a label for it
#         row = i + 1 # Start from the second row in the controls grid
#         toggle = Toggle(controls[row, 1], active=true)
#         Label(controls[row, 2], "Label $(label)", halign=:left)

#         # Link the toggle's state to the plot's visibility
#         # on(...) creates a listener that triggers when the toggle's state changes
#         on(toggle.active) do is_active
#             if haskey(point_cloud_plots, label)
#                 point_cloud_plots[label].visible = is_active
#             end
#         end
#     end
    
#     # Adjust column widths for a tidy layout
#     colsize!(f.layout, 2, Fixed(200)) # Give the controls a fixed width
#     rowgap!(controls, 10) # Add some space between the toggles

#     return f
# end

# function visualize_point_cloud_and_interface(
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int},
#     radii::Vector{Float64};
#     point_cloud_colormap = :Accent_6,
#     interface_colormap = :viridis,
#     show_wireframe = false,
# )
#     f = Figure(fontsize=12)

#     ax_mol = LScene(f[1, 1], show_axis=false)#, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))
#     ax_int = LScene(f[1, 2], show_axis=false, scenekw=(lights=[AmbientLight(RGBf(1.0, 1.0, 1.0))],))

#     _populate_point_cloud_scene!(
#         ax_mol, points, labels, radii,
#         colormap=point_cloud_colormap
#     )

#     interface_plot = _populate_interface_scene!(
#         ax_int, points, labels, radii,
#         show_wireframe=show_wireframe,
#         interface_colormap=interface_colormap
#     )

#     #cameracontrols!(ax_int.scene, cameracontrols(ax_mol.scene))
#     return f
# end


# #================================================================================#
# #                         INTERNAL HELPER FUNCTIONS                              #
# #================================================================================#

# function _compute_and_prepare_mesh_data(points::Vector{Vector{Float64}}, labels::Vector{Int}, radii::Vector{Float64}; center_point=nothing)
#     bcs, filtration = get_barycentric_subdivision_and_filtration(points, labels, radii)
    
#     # Extract faces (triangles) and vertex colors from the filtration data
#     faces = [TriangleFace(e[1]) for e in filtration if length(e[1]) == 3]
#     mesh_colors = [e[2] for e in filtration if length(e[1]) == 1]
    
#     # Convert barycentric coordinates to Point3f
#     bcs_points = [Point3f(bc) for bc in bcs]
    
#     # Center the points if a center_point is provided
#     if !isnothing(center_point)
#         bcs_points .-= center_point
#     end
    
#     mesh = GeometryBasics.Mesh(bcs_points, faces)
    
#     return mesh, mesh_colors, bcs_points
# end

# function _populate_point_cloud_scene!(
#     lscene::LScene, 
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int},
#     radii::Vector{Float64}; 
#     point_cloud_colormap = :Accent_6
# )
#     if isempty(points) || !(length(points) == length(radii) == length(labels))
#         return
#     end

#     meshscatter!(lscene, 
#         [Point3f(p) for p in points], 
#         marker = Sphere(Point3f(0), 1),
#         markersize = radii,
#         color = labels,
#         colormap = point_cloud_colormap
#     )
# end

# function _populate_interface_scene!(
#     lscene::LScene, 
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int}, 
#     radii::Vector{Float64};
#     show_wireframe = false, 
#     show_barycenters = false,
#     show_mc_points = false,
#     show_mc_edges = false,
#     interface_colormap = :viridis
# )
#     # 1. Compute mesh data (no centering for this static view)
#     msh, clr, bcs_points = _compute_and_prepare_mesh_data(points, labels, radii)

#     # 2. Determine color range for the mesh
#     min_v, max_v = isempty(clr) ? (0.0, 1.0) : (minimum(clr), maximum(clr))
    
#     # 4. Plot mesh and CAPTURE the plot object for the colorbar
#     mesh_plot = mesh!(lscene, msh, color=clr, colorrange=(min_v, max_v), colormap=interface_colormap, shading = false)
#     if show_wireframe
#         wireframe!(lscene, msh, color=:white, linewidth=1)
#     end

#     if show_barycenters
#         _populate_scene_with_barycenters!(lscene, bcs_points)
#     end

#     if show_mc_edges || show_mc_points
#         _populate_scene_with_multicolored_data(lscene, points, labels, radii; show_mc_edges=show_mc_edges, show_mc_points=show_mc_points)
#     end
# end

# function _populate_scene_with_barycenters!(scene, bcs_points; markersize=15)
#     scatter!(scene, bcs_points, color=:black, markersize=markersize, overdraw=true)
#     # Add text labels for each barycenter index
#     for (i, bc) in enumerate(bcs_points)
#         text!(scene, bc; text="$(i)", fontsize=15, color=:red, overdraw=true, align=(:center, :center))
#     end
# end

# function _populate_scene_with_multicolored_data(
#     lscene::LScene, 
#     points::Vector{Vector{Float64}}, 
#     labels::Vector{Int}, 
#     radii::Vector{Float64};
#     show_mc_points::Bool = false,
#     show_mc_edges::Bool = false
# )
#     color_map = let
#         cm = cgrad(:Accent_6, 6, categorical=true)
#         Dict(1 => cm[1], 2 => cm[2], 3 => cm[3], 4 => cm[6])
#     end
    
#     if show_mc_points
#         # Plot multicolored points
#         scatter!(lscene, [Point3f(p) for p in points], color=[color_map[l] for l in labels], strokecolor=:black, strokewidth=1)
#     end

#     if show_mc_edges
#         mc_tets = length(radii) == length(labels) ? get_multichromatic_tetrahedra(points, labels, radii) : get_multichromatic_tetrahedra(points, labels)
#         edge_points, edge_colors = _prepare_mc_edge_visuals(points, labels, mc_tets, color_map)
        
#         linesegments!(lscene, edge_points, color=edge_colors, linewidth=2)
#     end
# end
