using Random 
using GLMakie

function display_random_point_cloud_interface(number_of_points::Int = 10)
    points = [Vector{Float64}(e) for e in eachcol(rand(3, number_of_points))]
    color_labels = [rand(1:4) for _ in points]
    radii = [0.5 for _ in points]
    interface_surface = TopologicalInterfaces.InterfaceSurface(points, color_labels, radii)
    display(TopologicalInterfaces.get_interface_figure(interface_surface, points, color_labels, radii; show_wireframe = true, show_multicolored_points = true, show_multicolored_edges = true))
end