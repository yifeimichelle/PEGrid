using Base.Test
include("src/framework.jl")
framework = Framework("IRMOF-1")
@test_approx_eq_eps framework.crystaldensity() 698.1128626806609 .001
framework.chemicalformula()
@printf("framework constructed\n\n")

include("src/forcefield.jl")
forcefield = Forcefield("UFF", "Xe")
@printf("forcefield constructed\n\n")

include("src/energyutils.jl")

@test_approx_eq_eps energy_of_adsorbate("CH4", [.5,.2,.3], "IRMOF-1", "UFF", cutoff=12.5) [412049.736961] .1
fractional_coords = [.2 .3; .1 .4; .5 .8]
E = energy_of_adsorbate("CH4", fractional_coords, "IRMOF-1", "UFF", cutoff=12.5)
E_CH2CH2 = energy_of_adsorbate("CH2CH2", fractional_coords, "IRMOF-1", "UFF", cutoff=12.5)

# energy grid
@printf("energy grid test...\n\n\n")
include("src/energygrid.jl")
writegrid("CH4", "IRMOF-1", "UFF", gridspacing=1.0, cutoff=12.5)

# surface area
@printf("surface area test\n")
include("src/surface_area.jl")
sa = surface_area("N2", "IRMOF-1", "UFF", probe_size=3.31, num_sampling=10000)
