
using PDBTools

root = @__DIR__

function get_protor_radii(pdb_file::String, chain::String)
    py"""
    import freesasa
    import numpy as np
    def get_protor_radii(pdb_file, chain):
        print(pdb_file)
        structure = freesasa.Structure(pdb_file)
        return np.array([structure.radius(i) for i in range(structure.nAtoms()) if structure.chainLabel(i) == chain]) 
    """
    return py"get_protor_radii"(pdb_file, chain)
end

function get_template_centers_and_radii(mol_type::String)
    return default_generation(mol_type)
end

function default_generation(mol_type::String; chn = "A")
    file = "$(root)/pdbs/$(mol_type).pdb"
    atoms_A = select(read_pdb(file, "chain $(chn)"), a -> element(a) != "H" && resname(a) != "HOH")
    moveto!(atoms_A)
    template_centers = reduce(hcat,[Vector{Float64}(e) for e in coor(atoms_A)])
    template_radii = get_protor_radii(file, chn)
    template_centers, template_radii
end