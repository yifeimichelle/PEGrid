# push!(LOAD_PATH, homedir() * "/Dropbox/PEGrid/src")
module PEGrid

# this is where crystal structure files, adsorbate and forcefield input files are stored
global const PEGRID_DATA_DIR = homedir() * "/Dropbox/PEGrid/data"
 # global const PEGRID_OUTPUT_DIR = homedir() * "/PEGrid_output"

include("henry.jl")
include("gcmc.jl")
include("energygrid.jl")
include("getEwaldparams.jl")
export Framework, # framework.jl
       Adsorbate, uniform_unit_vector_on_sphere, random_rotation_matrix, # adsorbate.jl
       Forcefield, # forcefield.jl
       henry, # henry.jl
       # energyutils.jl
       electrostatic_potential, vdW_energy_of_bead, vdW_energy_of_adsorbate, get_replication_factors, electrostatic_energy_of_adsorbate, _generate_epsilons_sigmas2,
       # energygrid.jl
       write_vdW_grid, write_electrostatic_grid, # Grid conflicts with Gadfly...
       getEwaldparams, # getEwaldparams.jl
       grand_canonical_mc, adsorption_isotherm # gcmc.jl
end
