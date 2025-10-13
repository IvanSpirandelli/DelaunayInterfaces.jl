using Random 
using GLMakie

function display_random_point_cloud_interface(number_of_points::Int = 10, box_length = 1.0, radius = 0.5; show_wireframe = true, show_multicolored_points = true, show_multicolored_edges = true)
    points = [Vector{Float64}(e) for e in eachcol(rand(3, number_of_points).*box_length)]
    color_labels = [rand(1:4) for _ in points]
    radii = [radius for _ in points]
    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels, radii)
    display(TopologicalInterfaces.get_interface_figure(interface_surface, points, color_labels, radii; 
    show_wireframe = show_wireframe, show_multicolored_points = show_multicolored_points, show_multicolored_edges = show_multicolored_edges))
end