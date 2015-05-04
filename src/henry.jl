###
#   Compute Henry constant
#   Author: Cory M Simon CoryMSimon@gmail.com
###
include("energyutils.jl")

function henry(adsorbatename::String, structurename::String, forcefieldname::String, temperature::Float64; insertions_per_A3::Int=750, cutoff::Float64=12.5)
    """
    Compute Henry constant and ensemble average energy of an adsorbate inside a structure via Widom insertions
    """
    @printf("Constructing framework object for %s...\n", structurename)
    framework = Framework(structurename)
    
    @printf("Constructing adsorbate %s...\n", adsorbatename)
    adsorbate = Adsorbate(adsorbatename)
    if (adsorbate.nbeads > 1) & (temperature == -1.0)
        error("Provide temperature for Boltzmann weighted rotations, nbeads > 1 in adsorbate.")
    end

    @printf("Constructing forcefield(s) for bead(s) in %s...\n", forcefieldname)
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        @printf("\tBead %s...\n", adsorbate.bead_names[b])
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get array of framework atom positions and corresponding epsilons and sigmas for speed
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    @printf("Unit cell replication factors for LJ cutoff of %.2f A: %d by %d by %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    
    num_insertions = convert(Int, ceil(framework.v_unitcell * insertions_per_A3))
    @printf("Performing %d Widom insertions...\n", num_insertions)

    boltzmann_factor_sum = 0.0
    boltzmann_weighted_energy_sum = 0.0
    for i = 1:num_insertions
        xf_insert = [rand(), rand(), rand()]

        # translate adsorbate to this point
        adsorbate.translate(framework.f_to_cartesian_mtrx * xf_insert)
        # select a random orientiation if applicable
        if adsorbate.nbeads > 1
            adsorbate.perform_uniform_random_rotation()
        end
        
        # compute energy here
        _energy = _energy_of_adsorbate!(adsorbate,
                                    epsilons,
                                    sigmas,
                                    framework,
                                    rep_factors, 
                                    cutoff)
        boltzmann_weight = exp(-_energy / temperature)
        boltzmann_factor_sum += boltzmann_weight
        boltzmann_weighted_energy_sum += boltzmann_weight * _energy

        # translate adsorbate back to COM
        adsorbate.set_origin_at_COM()
    end
    
    KH = boltzmann_factor_sum / num_insertions / 8.314 / temperature
    avg_E = boltzmann_weighted_energy_sum / boltzmann_factor_sum
    @printf("%s Henry constant in %s at %.1f K = %f (mol/(m3-Pa))\n", adsorbatename, structurename, temperature, KH)
    @printf("\t<E> = %f K = %f kJ/mol\n", avg_E, avg_E * 8.314 / 1000.0)
end
