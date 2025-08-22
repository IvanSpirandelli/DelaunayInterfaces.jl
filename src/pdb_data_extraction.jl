function get_points_radii_labels(pdb_file::String, chain_ids::Vector{String})
    # 1. Get all atom properties, including atom names
    protein_data = get_atom_properties(pdb_file, chain_ids)

    # 2. Initialize the vectors to hold data from all chains
    all_points = Vector{Vector{Float64}}()
    all_radii = Float64[]
    all_labels = Int[]
    all_atom_names = String[] # <-- We now need to collect atom names

    # 3. Aggregate data from the requested chains
    for (label, chain_id) in enumerate(chain_ids)
        if haskey(protein_data, chain_id)
            chain_data = protein_data[chain_id]
            
            append!(all_points, chain_data.coords)
            append!(all_radii, chain_data.radii)
            append!(all_atom_names, chain_data.atom_names) # <-- Collect the names
            
            n_atoms_in_chain = length(chain_data.radii)
            append!(all_labels, fill(label, n_atoms_in_chain))
        end
    end

    if isempty(all_points)
        return all_points, all_radii, all_labels
    end

    # 4. Calculate the true center of mass using our new function
    center_of_mass = calculate_center_of_mass(all_points, all_atom_names)

    # 5. Subtract the center of mass from every point
    centered_points = [p .- center_of_mass for p in all_points]
    
    # 6. Return the centered points and other data
    return centered_points, all_radii, all_labels
end

function get_atom_properties(pdb_file::String, chain_ids::Vector{String})
    structure = fs[].Structure(pdb_file)
    results = Dict{String, NamedTuple}()
    n_atoms = pyconvert(Int, structure.nAtoms())

    for chain_id in chain_ids
        results[chain_id] = (
            coords = Vector{Vector{Float64}}(), # <-- ADDED
            radii = Float64[],
            atom_names = String[],
            residue_names = String[],
            residue_numbers = Int[]
        )
    end

    for i in 0:(n_atoms - 1)
        current_chain = pyconvert(String, structure.chainLabel(i))

        if haskey(results, current_chain)
            push!(results[current_chain].coords, pyconvert(Vector{Float64}, structure.coord(i)))
            push!(results[current_chain].radii, pyconvert(Float64, structure.radius(i)))
            push!(results[current_chain].atom_names, pyconvert(String, structure.atomName(i)))
            push!(results[current_chain].residue_names, pyconvert(String, structure.residueName(i)))
            push!(results[current_chain].residue_numbers, parse(Int, pyconvert(String, structure.residueNumber(i))))
        end
    end

    return results
end


const ATOM_MASSES = Dict(
    "H"  => 1.008,
    "C"  => 12.011,
    "N"  => 14.007,
    "O"  => 16.000,
    "S"  => 32.065,
    "P"  => 30.974,
)

"""
    get_element_symbol(atom_name::String) -> String

Extracts the chemical element symbol from a PDB atom name (e.g., "CA" -> "C").
This simple version assumes the element is the first non-numeric character.
"""
function get_element_symbol(atom_name::String)
    # Find the first character that is a letter and return it as a string
    for char in atom_name
        if isletter(char)
            return uppercase(string(char))
        end
    end
    return "C" # Fallback to Carbon if no letter is found
end

function calculate_center_of_mass(points::Vector{Vector{Float64}}, atom_names::Vector{String})
    if length(points) != length(atom_names)
        error("The number of points must match the number of atom names.")
    end

    if isempty(points)
        return [0.0, 0.0, 0.0] 
    end

    total_mass = 0.0
    weighted_coord_sum = [0.0, 0.0, 0.0]

    for (i, p) in enumerate(points)
        atom_name = atom_names[i]
        element = get_element_symbol(atom_name)
        mass = get(ATOM_MASSES, element, 1.0)
        total_mass += mass
        weighted_coord_sum .+= p .* mass
    end
    if total_mass == 0.0
        return [0.0, 0.0, 0.0]
    end
    return weighted_coord_sum / total_mass
end