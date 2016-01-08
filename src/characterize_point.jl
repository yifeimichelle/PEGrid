include("energyutils.jl")


function characterize_point(structurename::AbstractString, 
                        forcefieldname::AbstractString, 
                        adsorbate::Adsorbate;
                        r_bins::Array{Float64}=linspace(0, 12.5),
                        cutoff::Float64=12.5,
                        verbose_flag::Bool=false)
    """
    For a given adsorbate molecule position:
    -- find the contribution to the energy of each atom type within the cutoff
    -- get radial distribution function of each atom in the form of counts of atoms within radial bins (r_bins)
    -- get list of atoms, distances from binding site, and energies
    """
    if (verbose_flag)
        @printf("Constructing framework object for %s...\n", structurename)
    end
    framework = Framework(structurename)

    if (verbose_flag)
        @printf("Constructing forcefield(s) for bead(s) in %s...\n", forcefieldname)
    end
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        if (verbose_flag)
            @printf("\tBead %s...\n", adsorbate.bead_names[b])
        end
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework, cutoff)
    @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    
    # get position array and epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)

    # get unique atoms in framework
    uniqueatoms = unique(framework.atoms)
    println(uniqueatoms)
    
    # create dictionaries that keep track of energy contributions and rdf by each atom in framework
    energy_contributions = Dict() 
    rdf_by_atom = Dict()
    for a = 1:length(uniqueatoms)
        energy_contributions[uniqueatoms[a]] = 0.0
        rdf_by_atom[uniqueatoms[a]] = zeros(Int, length(r_bins) - 1)
    end

    # file for atoms, their distance from binding site, and energies
    if ! isdir(homedir() * "/PEGrid_output/bindingsites")
       mkdir(homedir() * "/PEGrid_output/bindingsites")
    end
    surrounding_atom_file = open(homedir() * "/PEGrid_output/bindingsites/" * structurename * "_surrounding_atom_list.csv", "w")
    @printf(surrounding_atom_file, "Atom,r(A),E(kJ/mol)\n")
   
    # Loop over all interactions
    E = 0.0  # initialize energy at this grid point
    for b = 1:adsorbate.nbeads  # for each bead
        # get fractional coordinate of bead
        fractional_coord = framework.cartesian_to_f_mtrx * adsorbate.bead_xyz[:, b]
        for i = 1:3 
            fractional_coord[i] = mod(fractional_coord[i], 1.0)
        end
        
        for rep_x = -rep_factors[1]:rep_factors[1]
            for rep_y = -rep_factors[2]:rep_factors[2]
                for rep_z = -rep_factors[3]:rep_factors[3]
                    for i_a = 1:framework.natoms
                        dx = framework.fractional_coords[:, i_a] + 1.0 * [rep_x, rep_y, rep_z] - fractional_coord
                        dx = framework.f_to_cartesian_mtrx * dx
                        r2 = dot(dx, dx)
                        if r2 < cutoff * cutoff
                            # bin this atom's radius in r_bins for radial distn function
                            bin_number = searchsortedfirst(r_bins, sqrt(r2)) - 1
                            if (r2 < maximum(r_bins) ^ 2) & (r2 != 0.0)
                                rdf_by_atom[framework.atoms[i_a]][bin_number] += 1
                            end

                            # compute VdW energy with Lennard-Jones potential.
                            sig_ovr_r6 = sigmas[i_a] * sigmas[i_a] ./ r2  # (sigma / r )^2
                            sig_ovr_r6 = sig_ovr_r6 * sig_ovr_r6 * sig_ovr_r6
                            E_here = 4.0 * epsilons[i_a] .* sig_ovr_r6 .* (sig_ovr_r6 - 1.0)
                            
                            # keep track of which atom this is from
                            energy_contributions[framework.atoms[i_a]] += E_here
                            
                            # write to rdf by atom file
                            @printf(surrounding_atom_file, "%s,%f,%f\n", framework.atoms[i_a], sqrt(r2), E_here * 8.314/ 1000.0)

                            # update total energy
                            E += E_here
                        end  # end if r2 within cutoff
                    end  # end loop over framework atoms
                end  # end replication in x-direction
            end  # end replication in y-direction
        end  # end replication in z-direction
    end  # loop over beads
    
    # assert total energy is accounted for
    E_total_check = 0.0
    for a = 1:length(uniqueatoms)
        E_total_check += energy_contributions[uniqueatoms[a]]
    end
    @assert(E > E_total_check - 0.00001, @sprintf("E=%f, E total in contributions = %f\n", E, E_total_check))
    @assert(E < E_total_check + 0.00001, @sprintf("E=%f, E total in contributions = %f\n", E, E_total_check))

    # write Energy contribution results to file
    if ! isdir(homedir() * "/PEGrid_output/energycontributions")
       mkdir(homedir() * "/PEGrid_output/energycontributions") 
    end
    energyfile = open(homedir() * "/PEGrid_output/energycontributions/" * structurename * "_" * adsorbate.name * "_" * forcefieldname * ".csv", "w")
    @printf(energyfile, "Binding site at fractional point: [ ")
    for b = 1:adsorbate.nbeads
        @printf(energyfile, "(%f, %f, %f) ", adsorbate.bead_xyz[1, b],  adsorbate.bead_xyz[2, b], adsorbate.bead_xyz[3, b])
    end
    @printf(energyfile, "]\n")
    @printf(energyfile, "Total energy = %f K\n", E)
    @printf(energyfile, "Forcefield: %s\n", forcefieldname)
    @printf(energyfile, "Adsorbate: %s\n", adsorbate.name)
    @printf(energyfile, "Atom,EnergyContribution(K)\n")

    for a = 1:length(uniqueatoms)
        @printf(energyfile, "%s,%f\n", uniqueatoms[a], 
            energy_contributions[uniqueatoms[a]])
    end
    close(energyfile)
   
    # write radial distn function results to file
    if ! isdir(homedir() * "/PEGrid_output/RDF")
       mkdir(homedir() * "/PEGrid_output/RDF") 
    end
    rdffile = open(homedir() * "/PEGrid_output/RDF/" * structurename * ".csv", "w")
    @printf(rdffile, "Binding site at fractional point: [ ")
    for b = 1:adsorbate.nbeads
        @printf(rdffile, "(%f, %f, %f) ", adsorbate.bead_xyz[1, b],  adsorbate.bead_xyz[2, b], adsorbate.bead_xyz[3, b])
    end
    @printf(rdffile, "]\n")
    
    # header
    @printf(rdffile, "bin_end")
    for a = 1:length(uniqueatoms)
        @printf(rdffile, ",%s", uniqueatoms[a])
    end
    @printf(rdffile, "\n")
    for b = 1:length(r_bins)-1
        @printf(rdffile, "%f", r_bins[b+1])
        for a = 1:length(uniqueatoms)
            @printf(rdffile, ",%f", rdf_by_atom[uniqueatoms[a]][b])
        end
        @printf(rdffile, "\n")
    end
    close(rdffile)
    
    return energy_contributions, rdf_by_atom
end

function find_radius_where_most_energy_is(x_f::Array{Float64},
                                          fraction::Float64,
                                          structurename::AbstractString, 
                                          forcefieldname::AbstractString, 
                                          adsorbate::Adsorbate;
                                          cutoff::Float64=12.5,
                                          verbose_flag::Bool=false)
    """
    How far from point x_f do we need to go to capture most of the energy? cutoff.
    """
    @assert((fraction > 0.0) & (fraction < 1.0))

    framework = Framework(structurename)

    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        if (verbose_flag)
            @printf("\tBead %s...\n", adsorbate.bead_names[b])
        end
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework, cutoff)
    
    # get epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)

    # get energy at this point
    E = _vdW_energy_of_adsorbate!(adsorbate,
            epsilons, 
            sigmas,
            framework,
            rep_factors,
            cutoff)
    # define array of r to search
    cutoffs = linspace(1.0, cutoff, 150)

    for i = 1:length(cutoffs)
        E_here = _vdW_energy_of_adsorbate!(adsorbate,
                epsilons, 
                sigmas,
                framework,
                rep_factors,
                cutoffs[i])
        if abs(E_here) >= abs(E) * fraction
            @printf("\tE at cutoff %f A = %f\n", cutoffs[i], E_here*8.314/1000.0)
            @printf("\tE at cutoff %f A = %f\n", cutoff, E*8.314/1000.0)
            return cutoffs[i]
        end
    end
end

function write_xyz_of_sphere(structurename::AbstractString, 
            x_f::Float64, y_f::Float64, z_f::Float64, 
            r::Float64; adsorbate=None)
    """
    Write xyz file of atoms in framework about a sphere of radius r 
    at fractional coords (x_f, y_f, z_f)
    pass name of adsorbate if want adsorbate position at (x_f, y_f, z_f) printed as well
    """
    if ! isdir(homedir() * "/PEGrid_output/bindingsites")
       mkdir(homedir() * "/PEGrid_output/bindingsites") 
    end
    if (adsorbate == None)
        xyz_file = open(homedir() * "/PEGrid_output/bindingsites/" * structurename * ".xyz", "w")
    else
        xyz_file = open(homedir() * "/PEGrid_output/bindingsites/" * structurename * "_" * adsorbate * ".xyz", "w")
    end

    framework = Framework(structurename)

    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework, r)
  
    # store atoms inside the sphere here
    atomtype = AbstractString[]
    x_ = Float64[]
    y_ = Float64[]
    z_ = Float64[]

    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3] 
                
                for a = 1:framework.natoms
                    x_atom = framework.fractional_coords[:, a] + 1.0 * [rep_x, rep_y, rep_z]  # fractional coord of atom
                    # vector between site (x_f, y_f, z_f) and atom in this unit cell
                    dx = [x_f, y_f, z_f] - x_atom

                    # convert vector to Cartesian coords
                    dx = framework.f_to_cartesian_mtrx * dx
                
                    # compute distance squared between grid point and each framework atom, r2
                    r2 = sum(dx .* dx)

                    # if within sphere, add to list to print to xyz
                    if r2 < r * r
                        # convert vector to Cartesian coords
                        x_atom = framework.f_to_cartesian_mtrx * x_atom
                        # add to list
                        push!(x_, x_atom[1])
                        push!(y_, x_atom[2])
                        push!(z_, x_atom[3])
                        push!(atomtype, framework.atoms[a]) 
                    end
                end  # end loop over framework atoms
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction

    # write atoms within sphere to xyz file
    n = length(x_)  # count of atoms in this sphere
    if adsorbate == None
        @printf(xyz_file, "%d\n\n", n)
    else  
        # write adsorbate position at x_f, y_f, z_f
        x = framework.f_to_cartesian_mtrx * [x_f, y_f, z_f]
        @printf(xyz_file, "%d\n\n", n + 1)
        @printf(xyz_file, "%s %f %f %f\n", adsorbate, x[1], x[2], x[3])
    end

    for i = 1:n
        @printf(xyz_file, "%s %f %f %f\n", atomtype[i],
                x_[i], y_[i], z_[i])
    end
    
    close(xyz_file)
end
