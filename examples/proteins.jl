using GLMakie
using JLD2

function get_example_data(example_id = "4bmg_dimer")
    @load "examples/data/$(example_id).jld2" points color_labels radii
    points, color_labels, radii
end

function display_protein_protein_interface(example_id = "4bmg_dimer"; alpha = true, solvent_radius = 1.4, show_wireframe = false)
    points, color_labels, radii = get_example_data(example_id)
    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels, radii .+ solvent_radius, alpha)
    display(TopologicalInterfaces.get_interface_figure(interface_surface; show_wireframe = show_wireframe))
end

function display_protein_space_filling_diagram(example_id = "4bmg_dimer")
    points, color_labels, radii = get_example_data(example_id)
    display(TopologicalInterfaces.get_point_cloud_figure(points, color_labels, radii; point_cloud_colormap = :magma))
end

function display_protein_protein_interface_and_space_filling_diagram(example_id = "4bmg_dimer"; alpha = true, solvent_radius = 1.4, show_wireframe = false)
    points, color_labels, radii = get_example_data(example_id)
    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels, radii .+ solvent_radius, alpha)
    display(TopologicalInterfaces.get_interface_and_point_cloud_figure(interface_surface, points, color_labels, radii; show_wireframe = show_wireframe, point_cloud_colormap = :magma))
end

function display_protein_protein_interface_filtration(example_id = "4bmg_dimer"; alpha = true, solvent_radius = 1.4, show_wireframe = false)
    points, color_labels, radii = get_example_data(example_id)
    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels, radii .+ solvent_radius, alpha)
    display(TopologicalInterfaces.get_interface_filtration_figure(interface_surface; show_wireframe = show_wireframe))
end
