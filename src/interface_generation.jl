using Distances

function get_multicolored_tetrahedra(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64}, weighted::Bool = true, alpha::Bool = true)
    if weighted
        if alpha
            get_multicolored_tetrahedra_of_weighted_alpha_complex(points, color_labels, radii)
        else
            get_multicolored_tetrahedra_of_weighted_delaunay_complex(points, color_labels, radii)
        end
    else
        get_multicolored_tetrahedra_of_delaunay_complex(points, color_labels)
    end
end

function get_multicolored_tetrahedra_of_delaunay_complex(points::Vector{Vector{Float64}}, color_labels::Vector{Int})
    py"""
    import numpy as np
    import diode
    def get_multicolored_tetrahedra_of_delaunay_complex(points, color_labels):
        points_np = np.asarray(points)

        # Use a lambda function for a concise check
        is_multi = lambda vs: len(set(color_labels[v] for v in vs)) >= 2

        # Use Python's 'and' operator and 'len()' function
        tetrahedra_0_based = [
            vs for (vs, fs) in diode.fill_alpha_shapes(points_np)
            if len(vs) == 4 and is_multi(vs)
        ]
        # Convert the final list of tetrahedra to 1-based indexing for use in Julia
        return [[i + 1 for i in t] for t in tetrahedra_0_based]
    """
    [e for e in eachrow(py"get_multicolored_tetrahedra_of_delaunay_complex"(hcat(points...)', color_labels))]
end

function get_multicolored_tetrahedra_of_weighted_delaunay_complex(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64})
    py"""
    import numpy as np
    import diode
    def get_multicolored_tetrahedra_of_weighted_delaunay_complex(points, color_labels):
        points_np = np.asarray(points)

        # Use a lambda function for a concise check
        is_multi = lambda vs: len(set(color_labels[v] for v in vs)) >= 2

        # Use Python's 'and' operator and 'len()' function
        tetrahedra_0_based = [
            vs for (vs, fs) in diode.fill_weighted_alpha_shapes(points_np)
            if len(vs) == 4 and is_multi(vs)
        ]
        
        return [[i + 1 for i in t] for t in tetrahedra_0_based]
    """
    weighted_points = [[p[1], p[2], p[3], r^2] for (p,r) in zip(points, radii)]
    [e for e in eachrow(py"get_multicolored_tetrahedra_of_weighted_delaunay_complex"(hcat(weighted_points...)', color_labels))]
end


function get_multicolored_tetrahedra_of_weighted_alpha_complex(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64})
    py"""
    import numpy as np
    import diode
    def get_multicolored_tetrahedra_of_weighted_alpha_complex(points, color_labels):
        points_np = np.asarray(points)

        # Use a lambda function for a concise check
        is_multi = lambda vs: len(set(color_labels[v] for v in vs)) >= 2

        # Use Python's 'and' operator and 'len()' function
        tetrahedra_0_based = [
            vs for (vs, fs) in diode.fill_weighted_alpha_shapes(points_np)
            if len(vs) == 4 and is_multi(vs) and fs <= 0
        ]
        
        return [[i + 1 for i in t] for t in tetrahedra_0_based]
    """
    weighted_points = [[p[1], p[2], p[3], r^2] for (p,r) in zip(points, radii)]
    [e for e in eachrow(py"get_multicolored_tetrahedra_of_weighted_alpha_complex"(hcat(weighted_points...)', color_labels))]
end

function get_chromatic_partitioning(tet, color_labels::Vector{Int})
    divs = [color_labels[v] for v in tet]
    parts = Dict(k => Vector{Int}([]) for k in unique(divs))
    for (i, d) in enumerate(divs)
        push!(parts[d], tet[i])
    end
    sort([v for v in values(parts)], by=length, rev = true)
end

function get_barycenter(points::Vector{Vector{Float64}}, vertices::Vector{Int}) 
    sum(points[vertices]) / length(vertices)
end

function get_or_create_bc_simplex_id_and_val!(partitioning::Vector{Vector{Int}}, uob_to_barycenter_simplices, points)
    simplex = Tuple(sort(vcat(partitioning...)))
    try
        (id, val) = uob_to_barycenter_simplices[simplex]
        return true, (id, val)
    catch
        id = length(uob_to_barycenter_simplices) + 1
        if length(partitioning) == 2
            part_one, part_two = partitioning
            val = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two)) #/ 2.0
        elseif length(partitioning) == 3
            part_one, part_two, part_three = partitioning
            a = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two))
            b = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_three))
            c = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_three))
            val = (a + b + c) / 3.0
        elseif length(partitioning) == 4
            part_one, part_two, part_three, part_four = partitioning
            a = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_two))
            b = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_three))
            c = euclidean(get_barycenter(points, part_one), get_barycenter(points, part_four))
            d = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_three))
            e = euclidean(get_barycenter(points, part_two), get_barycenter(points, part_four))
            f = euclidean(get_barycenter(points, part_three), get_barycenter(points, part_four))
            val = (a + b + c + d + e + f) / 6.0
        end
        uob_to_barycenter_simplices[simplex] = (id, val)
        return false, (id, val)
    end
end 

function _extend_barycenter_triangulation_scaffold!(barycenters, vertices, edges, points, mc_combinations, uob_to_barycenter_simplices)
    created = Vector{Bool}([])
    
    for uob_simplex in mc_combinations
        exists, (id, val) = get_or_create_bc_simplex_id_and_val!(uob_simplex, uob_to_barycenter_simplices, points)
        push!(vertices, (id, val))
        push!(created, !exists)
    end

    mc_combinations = [vcat(comb...) for comb in mc_combinations]

    uob_tri_indices = [i for (i, tri) in enumerate(mc_combinations) if length(tri) == 3]
    uob_tet_index = findfirst(x -> length(x) == 4, mc_combinations)


    uob_edge_indices = [i for (i, edge) in enumerate(mc_combinations) if length(edge) == 2]
    bc_of_edges = [get_barycenter(points, comb) for comb in mc_combinations[uob_edge_indices]]

    new_barycenters = [[0.0, 0.0, 0.0] for i in 1:length(mc_combinations)]
    new_barycenters[uob_edge_indices] = bc_of_edges

    bc_of_tris = [get_barycenter(new_barycenters, [i for i in uob_edge_indices if issubset(mc_combinations[i], comb)]) for comb in mc_combinations[uob_tri_indices]]
    new_barycenters[uob_tri_indices] = bc_of_tris
    
    bc_of_tet = get_barycenter(new_barycenters, [i for i in uob_edge_indices if issubset(mc_combinations[i], mc_combinations[uob_tet_index])])
    new_barycenters[uob_tet_index] = bc_of_tet

    for (i,bc) in enumerate(new_barycenters)
        if created[i]
            push!(barycenters, bc)
        end
    end

    for (i, j) in [(i,j) for (i, c1) in enumerate(mc_combinations) for (j, c2) in enumerate(mc_combinations) if c1 != c2 && issubset(c1, c2)]
        push!(edges, (sort!([vertices[i][1], vertices[j][1]]), minimum([vertices[i][2], vertices[j][2]])))
    end
end

function get_barycentric_subdivision_and_filtration(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64} = Vector{Float64}(), weighted::Bool = true, alpha::Bool = true)
    @assert length(color_labels) == length(points) "Each point must have a corresponding color_label"
    if weighted
        @assert length(radii) == length(color_labels) == length(points) "Each point must have an assigned radius to compute weighted complexes."
    end
    
    mc_tets = get_multicolored_tetrahedra(points::Vector{Vector{Float64}}, color_labels::Vector{Int}, radii::Vector{Float64}, weighted::Bool, alpha::Bool)
    barycenters = Vector{Vector{Float64}}([])
    filtration = Set{Tuple{Vector{Int32}, Float64}}([])
    uob_to_barycenter_simplices = Dict{Any, Any}()
    for vs in mc_tets
        parts = get_chromatic_partitioning(vs, color_labels)
        vertices = Vector{Tuple{Int, Float64}}([])
        edges = Vector{Tuple{Vector{Int}, Float64}}([])
        triangles = Vector{Tuple{Vector{Int}, Float64}}([])
        if length(parts) == 2
            part_one, part_two = parts
            if length(part_one) == length(part_two) == 2
                u, v = part_one
                x, y = part_two
                mc_combinations = [[[u],[x]], [[v],[x]], [[v],[y]], [[u],[y]], [[u,v],[x]], [[v], [x,y]], [[u,v], [y]], [[u], [x,y]], [[u,v],[x,y]]]
                _extend_barycenter_triangulation_scaffold!(barycenters, vertices, edges, points, mc_combinations, uob_to_barycenter_simplices)
                for (i, j, k) in [(9,1,5), (9,5,2), (9,2,6), (9,6,3), (9,3,7), (9,7,4), (9,4,8), (9,8,1)]
                    push!(triangles, ([vertices[i][1], vertices[j][1], vertices[k][1]], minimum([vertices[i][2], vertices[j][2], vertices[k][2]])))
                end
            elseif 3 == length(part_one) != length(part_two) == 1
                u,v,w = part_one[1], part_one[2], part_one[3]
                x = part_two[1]
                mc_combinations = [[[u],[x]], [[v],[x]], [[w],[x]], [[u,v],[x]], [[v,w], [x]], [[u,w],[x]], [[u,v,w],[x]]]                
                _extend_barycenter_triangulation_scaffold!(barycenters, vertices, edges, points, mc_combinations, uob_to_barycenter_simplices)
                for (i, j, k) in [(7,1,4), (7,4,2), (7,2,5), (7,5,3), (7,3,6), (7,6,1)]
                    push!(triangles, ([vertices[i][1], vertices[j][1], vertices[k][1]], minimum([vertices[i][2], vertices[j][2], vertices[k][2]])))
                end
            else
                error("Invalid partitioning: $vs -> $parts")
            end
        elseif length(parts) == 3
            part_one, part_two, part_three = parts
            a, b = part_one
            u = part_two[1]
            x = part_three[1]
            mc_combinations = [[[a],[u]], [[a],[x]], [[b],[x]], [[b],[u]], [[u],[x]], [[a],[u],[x]], [[a,b],[x]], [[b],[u],[x]], [[a,b],[u]], [[a,b],[u],[x]]]
            _extend_barycenter_triangulation_scaffold!(barycenters, vertices, edges, points, mc_combinations, uob_to_barycenter_simplices)
            for (i, j, k) in [(10,1,6), (10,6,5), (10,5,8), (10,8,4), (10,4,9), (10,9,1), (10,3,8), (10,6,2), (10,2,7), (10,7,3)]
                push!(triangles, ([vertices[i][1], vertices[j][1], vertices[k][1]], minimum([vertices[i][2], vertices[j][2], vertices[k][2]])))
            end
        elseif length(parts) == 4
            part_one, part_two, part_three, part_four = parts
            a = part_one[1]
            i = part_two[1]
            u = part_three[1]
            x = part_four[1]
            mc_combinations = [[[a],[i]], [[a],[u]], [[a],[x]], [[i],[u]], [[i],[x]], [[u],[x]], [[a],[i],[u]], [[a],[i],[x]], [[i],[u],[x]],  [[a],[u],[x]], [[a],[i],[u],[x]]]
            _extend_barycenter_triangulation_scaffold!(barycenters, vertices, edges, points, mc_combinations, uob_to_barycenter_simplices)
            for (i, j, k) in [(11,4,9), (11,9,5), (11,5,8), (11,8,1), (11,1,7), (11,7,4), (11,10,2), (11,6,10), (11,9,6), (11,2,7), (11,10,3), (11,3,8)]
                push!(triangles, ([vertices[i][1], vertices[j][1], vertices[k][1]], minimum([vertices[i][2], vertices[j][2], vertices[k][2]])))
            end
        end

        union!(filtration, Set([([v[1]], v[2]) for v in vertices]))
        union!(filtration, Set(edges))
        union!(filtration, Set(triangles))
    end
    barycenters, sort!(sort!(collect(filtration), by = x -> x[1]), by = x -> length(x[1]))
end