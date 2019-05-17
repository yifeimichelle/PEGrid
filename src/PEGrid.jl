# push!(LOAD_PATH, homedir() * "/Dropbox/PEGrid/src")
module PEGrid

# this is where crystal structure files, adsorbate and forcefield input files are stored
global const PEGRID_DATA_DIR = homedir() * "/LSMO/Michelle/Research/PEGrid/data"
@printf("Your data directory is set as:\n%s\n", PEGRID_DATA_DIR)
@printf("Change the variable PEGRID_DATA_DIR in PEGrid.jl if you wish to change this\n")
 # global const PEGRID_OUTPUT_DIR = homedir() * "/PEGrid_output"

include("henry.jl")
 # include("gcmc.jl")
include("energygrid.jl")
include("getEwaldparams.jl")
include("positiondistn.jl")
export Framework, # framework.jl
       Adsorbate, uniform_unit_vector_on_sphere, random_rotation_matrix, # adsorbate.jl
       Forcefield, # forcefield.jl
       henry, # henry.jl
       # energyutils.jl
       electrostatic_potential, vdW_energy_of_bead, vdW_energy_of_adsorbate, get_replication_factors, electrostatic_energy_of_adsorbate, _generate_epsilons_sigmas2,
       # energygrid.jl
       write_vdW_grid, write_electrostatic_grid, # Grid conflicts with Gadfly...
       getEwaldparams, # getEwaldparams.jl
       grand_canonical_mc, adsorption_isotherm, # gcmc.jl
       writeprobabilitydistncube, convert_pdb_to_xyz # positiondistn.jl
end
