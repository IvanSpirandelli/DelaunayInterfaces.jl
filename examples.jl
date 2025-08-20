# Standard library and package imports
using Pkg

# Activate the project environment to ensure all dependencies are met
Pkg.activate("Project.toml")
Pkg.instantiate()

# Import necessary libraries for the simulation
ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using TopologicalInterfaces
using GLMakie

#================================================================================#
#                         EXAMPLE MANAGEMENT FRAMEWORK                           #
#================================================================================#

"""
    Example

A struct to encapsulate all data and metadata related to a single visualization example.

# Fields
- `id::Symbol`: A unique identifier for the example (e.g., `:tetra_3_1`).
- `points::Vector{Vector{Float64}}`: The point cloud data.
- `labels::Vector{Int}`: The label for each point.
- `description::String`: An optional longer description of the example.
"""
struct Example
    id::Symbol
    points::Vector{Vector{Float64}}
    labels::Vector{Int}
    description::String
end

# --- Database of all available examples ---
# Using a constant Dict with a Symbol key provides fast, type-stable lookups.
const EXAMPLES = Dict{Symbol, Example}(
    :tetra_3_1 => Example(
        :tetra_3_1, 
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]],
        [1,1,1,2],
        "A simple tetrahedron with three points of one label and one of another."
    ),

    :tetra_2_2 => Example(
        :tetra_2_2,
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2]],
        [1,1,2,2],
        "The interface bisects the tetrahedron."
    ),

    :tetra_2_1_1 => Example(
        :tetra_2_1_1,
        [[0.1, 0.1, 0.1], [2.0, 0.0, 0.2], [1.0, 2.5, 0.0], [1.0, 1.4, 1.0]],
        [1,1,2,3],
        "An interface generated from three distinct labels."
    ),

    :tetra_1_1_1_1 => Example(
        :tetra_1_1_1_1,
        [[0.0, 0.0, 0.0], [2.0, 0.15, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 2.4]],
        [1,2,3,4],
        "A fully multicolored tetrahedron with four labels, showing the internal Voronoi-like structure."
    ),

    :two_colors_small => Example(
        :two_colors_small,
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.4], [1.0, 3.0, 1.2], [4.0, 2.0, 1.8], [2.5, 2.5, 2.5]],
        [1,1,1,2,2,2],
        "A small system with two groups of three points each."
    )
)

# --- Create a globally sorted list of example IDs ---
# This list is sorted first by the number of points, then by the number of unique labels.
const SORTED_EXAMPLE_IDS = sort(
    collect(keys(EXAMPLES)),
    by = id -> (length(EXAMPLES[id].points), length(unique(EXAMPLES[id].labels)))
)


#================================================================================#
#                             PUBLIC API FUNCTIONS                               #
#================================================================================#

"""
    list_examples()

Prints a formatted list of all available examples in a sorted order (simplest first).
"""
function list_examples()
    println("="^40)
    println("           Available Examples")
    println("="^40)
    # Iterate over the globally sorted list of IDs
    for id in SORTED_EXAMPLE_IDS
        ex = EXAMPLES[id]
        println("ID   : $(ex.id)")
        if !isempty(ex.description)
            println("Desc : $(ex.description)")
        end
        println("-"^40)
    end
end

"""
    show_example(id::Symbol; ...)

Displays a visualization of the specified example with a standard set of options.
"""
function show_example(
    id::Symbol;
    show_mc_edges::Bool = true,
    show_wireframe::Bool = true,
    show_barycenters::Bool = false,
    interface_colormap = :viridis
)
    # Check if the requested example exists.
    if !haskey(EXAMPLES, id)
        @error "Example with ID ':$id' not found. Please choose from the following:"
        list_examples()
        return
    end

    ex = EXAMPLES[id]
    
    # (FIXED) Print the example's ID, as the 'name' field was removed.
    println("Displaying example: :$(ex.id)")
    display(TopologicalInterfaces.get_interface_and_multicolored_tetrahedron_visualization(
        ex.points,
        ex.labels;
        show_mc_edges = show_mc_edges,
        show_wireframe = show_wireframe,
        show_barycenters = show_barycenters,
        interface_colormap = interface_colormap
    ))
end

"""
    example_flip_book(; example_ids = SORTED_EXAMPLE_IDS)

Creates an interactive sequence visualization ("flip book") of multiple examples,
sorted from simplest to most complex by default.
"""
function example_flip_book(; example_ids = SORTED_EXAMPLE_IDS)
    # Filter the examples based on the provided IDs
    examples_to_show = [EXAMPLES[id] for id in example_ids if haskey(EXAMPLES, id)]
    
    if isempty(examples_to_show)
        @error "No valid examples selected for the flip book."
        return
    end

    points_sequence = [ex.points for ex in examples_to_show]
    labels_sequence = [ex.labels for ex in examples_to_show]

    # Display the sequence visualization
    display(TopologicalInterfaces.visualize_interface_sequence(
        points_sequence,
        labels_sequence;
        show_wireframe = true,
        global_colorrange = false
    ))
end