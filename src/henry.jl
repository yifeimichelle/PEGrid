###
#   Compute Henry constant
#   Author: Cory M Simon CoryMSimon@gmail.com
###
include("energyutils.jl")

function henry(adsorbatename::String, 
               structurename::String, 
               forcefieldname::String, 
               temperature::Float64; 
               insertions_per_A3::Int=750, 
               cutoff::Float64=12.5, 
               verboseflag::Bool=false, 
               write_to_file=false)
    """
    Compute Henry constant and ensemble average energy of an adsorbate inside a structure via Widom insertions

    Keep track of lowest energy configuration
    """
    if verboseflag
        @printf("Constructing framework object for %s...\n", structurename)
    end
    framework = Framework(structurename)
    
    if verboseflag
        @printf("Constructing adsorbate %s...\n", adsorbatename)
    end
    adsorbate = Adsorbate(adsorbatename)
    adsorbate_min_config = deepcopy(adsorbate)  # preallocate
    if (adsorbate.nbeads > 1) & (temperature == -1.0)
        error("Provide temperature for Boltzmann weighted rotations, nbeads > 1 in adsorbate.")
    end
    
    if verboseflag
        @printf("Constructing forcefield(s) for bead(s) in %s...\n", forcefieldname)
    end
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        if verboseflag
            @printf("\tBead %s...\n", adsorbate.bead_names[b])
        end
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get array of framework atom positions and corresponding epsilons and sigmas for speed
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    if verboseflag
        @printf("Unit cell replication factors for LJ cutoff of %.2f A: %d by %d by %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    end
    
    num_insertions = convert(Int, ceil(framework.v_unitcell * insertions_per_A3))
    if verboseflag
        @printf("Performing %d Widom insertions...\n", num_insertions)
    end

    # for keeping track of minimum energy
    E_min = Inf  # pre-allocate E_min as inf

    boltzmann_factor_sum = 0.0  # \sum_i e^{-\beta E_i}
    boltzmann_weighted_energy_sum = 0.0  # \sum_i E_i e^{-\beta E_i}
    for i = 1:num_insertions
        xf_insert = [rand(), rand(), rand()]

        # translate adsorbate to this point
        adsorbate.translate_to(framework.f_to_cartesian_mtrx * xf_insert)
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
        
        # record for keeping track of avg's for K_H
        boltzmann_weight = exp(-_energy / temperature)
        boltzmann_factor_sum += boltzmann_weight 
        boltzmann_weighted_energy_sum += boltzmann_weight * _energy 

        # for calculating min energy pos
        if (_energy < E_min)
            E_min = _energy  # update minimum energy
            adsorbate_min_config = deepcopy(adsorbate)
        end
    end
    
    KH = boltzmann_factor_sum / num_insertions / 8.314 / temperature
    avg_E = boltzmann_weighted_energy_sum / boltzmann_factor_sum
    @printf("%s Henry constant in %s at %.1f K = %f (mol/(m3-Pa))\n", adsorbatename, structurename, temperature, KH)
    @printf("\t<E> %f K = %f kJ/mol\n", avg_E, avg_E * 8.314 / 1000.0)
    @printf("Minimum energy encountered = %f kJ/mol\n", E_min * 8.314 / 1000.0)

    # Write to file
    if write_to_file
        if ! isdir(homedir() * "/PEGrid_output/henries")
           mkdir(homedir() * "/PEGrid_output/henries") 
        end
        f = open(homedir() * "/PEGrid_output/henries/" * structurename * "_" * adsorbatename * "_henry.txt", "w")
        @printf(f, "T = %f K\n", temperature)
        @printf(f, "Framework density (kg/m3) = %f\n", framework.crystaldensity())
        @printf(f, "Insertions per A3: %d\n", insertions_per_A3)
        @printf(f, "<E> (kJ/mol) = %f\n", avg_E * 8.314 / 1000)
        @printf(f, "Min. E = %f kJ/mol\n", E_min * 8.314 / 1000.0)
        @printf(f, "KH, %s (mol/m3-Pa) = %e\n", adsorbatename, KH)
        close(f)
    end
    return KH, avg_E, adsorbate_min_config 
end
