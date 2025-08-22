module TopologicalInterfaces

using PythonCall

const diode = Ref{Py}()
const np = Ref{Py}()
const fs = Ref{Py}()

function __init__()
    diode[] = pyimport("diode")
    np[] = pyimport("numpy")
    fs[] = pyimport("freesasa")
end
include("interface_generation.jl")
include("interface_visualization.jl")
include("pdb_data_extraction.jl")
include("examples.jl")
end