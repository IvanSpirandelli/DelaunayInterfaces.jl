module DelaunayInterfaces

using PyCall

include("interface_generation.jl")

struct InterfaceSurface
    vertices::Vector{Vector{Float64}}
    filtration::Vector{Tuple{Vector{Int32}, Float64}}
    weighted::Bool
    alpha::Bool
end

function InterfaceSurface(points::Vector{Vector{Float64}}, color_labels::Vector{Int})
    vertices, filtration = get_barycentric_subdivision_and_filtration(points, color_labels, Vector{Float64}(), false, false)
    InterfaceSurface(vertices, filtration, false, false)
end

function InterfaceSurface(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64}, alpha::Bool = true)
    if length(radii) == 0
        return InterfaceSurface(points, color_labels)
    end
    vertices, filtration = get_barycentric_subdivision_and_filtration(points, color_labels, radii, true, alpha)
    InterfaceSurface(vertices, filtration, true, alpha)
end

include("interface_visualization.jl")

end