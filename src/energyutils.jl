###
#  Author: Cory M. Simon (CoryMSimon@gmail.com)
###
include("framework.jl")
include("forcefield.jl")
include("adsorbates.jl")


function _electrostatic_potential_of_bead(xyz_coord::Array{Float64},
                         framework::Framework,
                         rep_factors::Array{Int},
                         cutoff::Float64)
    """
    EWald sum to compute electrostatic potential at Cartesian point xyz_coord
    
    :param Framework framework: the structure
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: short-range cutoff distance
    """
    @assert(size(xyz_coord) == (3,))
    # get fractional coord of bead
    fractional_coord = framework.cartesian_to_f_mtrx * xyz_coord
    
    # reflect bead to [0,1] via periodic boundary conditions
    for i = 1:3
        fractional_coord[i] = mod(fractional_coord[i], 1.0)
    end

    # Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
    # 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
    # 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
    epsilon_0 = 4.7622424954949676e-7;  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)
    alpha = 0.025  # related to variance of Gaussians in Ewald sums

    ######
    ######  Short-range
    ######
    E_sr = 0.0  # short range energy
    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3]
                # loop over framework atoms
                for i = 1:framework.natoms
                    # get distance of this bead to framework atom
                    dx = fractional_coord - framework.fractional_coords[:, i] + 1.0 * [rep_x, rep_y, rep_z]
                    # transfer to Cartesian coords
                    dx = framework.f_to_cartesian_mtrx * dx
                    r = sqrt(dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3])

                    # add contribution to sr energy
                    if r < cutoff
                        E_sr += framework.charges[i] / (4.0 * pi * epsilon_0) / r * erfc(r * sqrt(alpha))
                    end
                end
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction
    @printf("Esr = %f\n", E_sr)

    ######
    ######  Long-range
    ######
    k_reps = [7, 7, 6]  # number of times to replicate in the k-space
    E_lr = 0.0
    for kx = -k_reps[1]:k_reps[1]
        for ky = -k_reps[2]:k_reps[2]
            for kz = -k_reps[3]:k_reps[3]
                if (kx == 0) & (ky == 0) & (kz == 0)
                    continue
                end
                # get the k vector
                k = kx * framework.reciprocal_lattice[:, 1] + ky * framework.reciprocal_lattice[:, 2] + kz * framework.reciprocal_lattice[:, 3]
                magnitude_k_2 = dot(k, k)
                
                E_lr_this_k = 0.0

                # loop ovr framework atoms
                for i = 1:framework.natoms
                    # Cartesian coord of framework atom
                    x_i = framework.f_to_cartesian_mtrx * framework.fractional_coords[:, i]
                    E_lr_this_k += framework.charges[i] * cos(dot(k, x_i - xyz_coord))
                end
                E_lr_this_k = E_lr_this_k / magnitude_k_2 * exp(- magnitude_k_2 / (4.0 * alpha))
                E_lr += E_lr_this_k
            end
        end
    end
    E_lr = E_lr / framework.v_unitcell / epsilon_0
    @printf("Elr = %f\n", E_lr)
    return E_lr + E_sr
end

function _vdW_energy_of_bead(xyz_coord::Array{Float64},
                         epsilons::Array{Float64},
                         sigmas::Array{Float64},
                         framework::Framework,
                         rep_factors::Array{Int},
                         cutoff::Float64)
    """
    Compute Van der Waals interaction energy via the LJ potential of a bead at xyz_coord (Cartesian!)
    
    :param Array{Float64} xyz_coord: Cartesian coordinate of bead. shape: (3,)
    :param Array{Float64} epsilons: Lennard-Jones epsilon parameters for bead with corresponding framework atoms. shape: (framework.natoms,)
    :param Array{Float64} sigmas: Lennard-Jones sigma parameters for bead with corresponding framework atoms. shape: (framework.natoms,)
    :param Framework framework: the structure
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: Lennard-Jones cutoff distance
    """
    @assert(size(xyz_coord) == (3,))
    # get fractional coord of bead
    fractional_coord = framework.cartesian_to_f_mtrx * xyz_coord

    # reflect bead to [0,1] via periodic boundary conditions
    for i = 1:3
        fractional_coord[i] = mod(fractional_coord[i], 1.0)
    end

    E = 0.0  # initialize energy of this bead as 0 to subsequently add pairwise contributions
    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3]
                # what follows is vectorized over the framework atoms in the primitive unit cell
                # subtract from each framework atom position the grid point
                dx = broadcast(-, framework.fractional_coords, fractional_coord + 1.0 * [rep_x, rep_y, rep_z])
                
                # transfer to Cartesian coords
                dx = framework.f_to_cartesian_mtrx * dx
                
                # compute distance squared between grid point and each framework atom, r2
                r2 = vec(sum(dx .* dx, 1))

                # select which interactions are within the cutoff distance
                idx_within_cutoff = r2 .< cutoff * cutoff

                # compute VdW energy with Lennard-Jones potential.
                sig_ovr_r6 = sigmas[idx_within_cutoff] .* sigmas[idx_within_cutoff] ./ r2[idx_within_cutoff]  # (sigma / r )^2
                sig_ovr_r6 = sig_ovr_r6 .* sig_ovr_r6 .* sig_ovr_r6
                
                E += sum(4.0 * epsilons[idx_within_cutoff] .* sig_ovr_r6 .* (sig_ovr_r6 - 1.0))
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction

    return E
end

function _vdW_energy_of_adsorbate!(adsorbate::Adsorbate,
            epsilons::Array{Float64}, 
            sigmas::Array{Float64},
            framework::Framework,
            rep_factors::Array{Int}, 
            cutoff::Float64)
    """
    Compute Van der Waals interaction energy via the Lennard Jones potential.
    Vectorized: this function takes an Array{Float64} of :math:`r^2` distances from each framework atom
    
    :param Adsorbate adsorbate: the molecule (includes position) that we want to calc the energy of
    :param Array{Float64} epsilons: Lennard-Jones epsilon parameters corresponding to each framework atom in the array
    :param Array{Float64} sigmas: Lennard-Jones sigma parameters corresponding to each framework atom in the array
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: Lennard-Jones cutoff distance
    :returns: Float64 E: energy of adsorbate at (x_f,y_f,z_f) in Kelvin (well, same units as epsilon)
    """
    E = 0.0
    for b = 1:adsorbate.nbeads  # for each bead
        E += _vdW_energy_of_bead(adsorbate.bead_xyz[:, b], epsilons[:, b], sigmas[:, b], framework, rep_factors, cutoff)
    end  # end loop over beads in adsorbate

    return E  # in Kelvin
end

function get_replication_factors(f_to_cartesian_mtrx::Array{Float64}, cutoff::Float64)
    """
    Get replication factors for the unit cell such that only directly adjacent ghost unit cells 
        are required for consideration in applying periodic boundary conditions (PBCs).
    That is, a particle in the "home" unit cell will not be within the cutoff distance of 
        any particles two ghost cells away from the home cell.
    This is done by ensuring that a sphere of radius cutoff can fit inside the unit cell
    
    :param Array{Float64} f_to_cartesian_mtrx: transformation matrix from fractional to Cartesian coordinates
    :param Float64 cutoff: Lennard-Jones cutoff distance, beyond which VdW interactions are approximated as zero.
    """
    # unit cell vectors
    a = f_to_cartesian_mtrx[:,1]
    b = f_to_cartesian_mtrx[:,2]
    c = f_to_cartesian_mtrx[:,3]

    # normal vectors to the faces of the unit cell
    n_ab = cross(a, b)
    n_bc = cross(b, c)
    n_ac = cross(a, c)

    # vector at center of unit cell. i.e center of sphere that we are ensuring fits within the unit cell
    c0 = [a b c] * [0.5, 0.5, 0.5]
    
    rep_factors = [1, 1, 1] # replication factors stored here. Initialize as 1 x 1 x 1 supercell.
    # replicate unit cell until the distance of c0 to the 3 faces is greater than cutoff/2
    # distance of plane to point: http://mathworld.wolfram.com/Point-PlaneDistance.html
    # a replication.
    while abs(dot(n_bc, c0)) / vecnorm(n_bc) < cutoff / 2.0
        rep_factors[1] += 1
        a = 2 * a
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # b replication
    while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff / 2.0
        rep_factors[2] += 1
        b = 2 * b
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # c replication
    while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff / 2.0
        rep_factors[3] += 1
        c = 2 * c
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    
    return rep_factors
end

function _generate_epsilons_sigmas(framework::Framework, forcefields::Array{Forcefield})
    """
    For speeding up energy computations, we use arrays of LJ parameters corresponding to framework atoms.

    Arrays `epsilons` and `sigmas` contain corresponding epsilon and sigma LJ parameters for each framework atom.
    """
    # array of LJ parameters that correspond to framework atoms
    n_beads = length(forcefields)
    epsilons = Array(Float64, (framework.natoms, n_beads)) # make same shape as r2...
    sigmas = Array(Float64, (framework.natoms, n_beads)) 
    for ff = 1:n_beads
        for i = 1:framework.natoms
            if ~ haskey(forcefields[ff].epsilon, framework.atoms[i])
                error(@sprintf("Atom %s not present in force field.", framework.atoms[i]))
            end
            epsilons[i, ff] = forcefields[ff].epsilon[framework.atoms[i]]
            sigmas[i, ff] = forcefields[ff].sigma[framework.atoms[i]]
        end
    end

    return  epsilons, sigmas
end

function vdW_energy_of_adsorbate_config(adsorbate::Adsorbate, structurename::String, forcefieldname::String; cutoff::Float64=12.5)
    """
    Get energy of adsorbate configuration (pass adsorbate object
    """
    framework = Framework(structurename)
    
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    
    # get position array and epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)
    E = _vdW_energy_of_adsorbate!(adsorbate,
                                epsilons,
                                sigmas,
                                framework,
                                rep_factors, 
                                cutoff)
    return E
end


function vdW_energy_of_adsorbate(adsorbatename::String,
                        fractional_translations::Array{Float64},
                        structurename::String,
                        forcefieldname::String;
                        cutoff::Float64=12.5,
                        num_rotation_samples::Int=500,
                        temperature::Float64=-1.0,
                        verbose_flag::Bool=false)
    """
    Compute energy of adsorbate

    :returns Array{Float64} of potential energies at each point in Kelvin (well, same units as \epsilon given in forcefield)
    """
    if (verbose_flag)
        @printf("Constructing adsorbate %s...\n", adsorbatename)
    end
    adsorbate = Adsorbate(adsorbatename)
    if (adsorbate.nbeads > 1) & (temperature == -1.0)
        error("Need to input temperature (in K) for adsorbate with >1 beads for Boltmann weighting of orientations (rotations)")
    end

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
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    if (verbose_flag)
        @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    end
    
    # get position array and epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)

    if (size(fractional_translations)[1] != 3)
        error("fractional_translations mtrx shld be 3 by npoints")
    end
    npoints = (length(fractional_translations) != 3) ? size(fractional_translations)[2] : 1

    if (verbose_flag)
        @printf("\tComputing potential energy at %d points\n", npoints)
    end

    E = zeros(npoints)  # pre-allocate array of energies

    for i = 1:npoints
        # translate adsorbate by vector in cartesian space:
        adsorbate.translate_to(framework.f_to_cartesian_mtrx * fractional_translations[:, i])

        # compute energy
        if adsorbate.nbeads == 1  # no need for sampling rotations
            E[i] = _vdW_energy_of_adsorbate!(adsorbate,
                                        epsilons,
                                        sigmas,
                                        framework,
                                        rep_factors, 
                                        cutoff)
        else  # need to sample rotations and get boltzmann-weighted average
            boltzmann_weight_sum = 0.0
            weighted_energy_sum = 0.0
            for k = 1:num_rotation_samples
                adsorbate.perform_uniform_random_rotation()
                _energy = _vdW_energy_of_adsorbate!(adsorbate,
                                            epsilons,
                                            sigmas,
                                            framework,
                                            rep_factors, 
                                            cutoff)
                boltzmann_weight = exp(-_energy / temperature)
                boltzmann_weight_sum += boltzmann_weight
                weighted_energy_sum += boltzmann_weight * _energy
            end
            E[i] = weighted_energy_sum / boltzmann_weight_sum
        end
    end
    
    return E  # in (Kelvin)
end

function find_min_energy_position(structurename::String, 
                                  forcefieldname::String, 
                                  adsorbatename::String,
                                  x_f_start::Array{Float64};
                                  num_rotation_samples::Int=200,
                                  temperature::Float64=-1.0,
                                  cutoff::Float64=12.5)
    """
    Find minimum energy position of adsorbate in framework
    """
    @printf("Constructing framework object for %s...\n", structurename)
    framework = Framework(structurename)

    @printf("Constructing adsorbate %s...\n", adsorbatename)
    adsorbate = Adsorbate(adsorbatename)
    if (adsorbate.nbeads > 1) & (temperature == -1.0)
        error("Need to provide temperature to calculate energy with different rotations\n")
    end
    
    @printf("Constructing forcefield object for %s...\n", forcefieldname)
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        @printf("\tBead %s...\n", adsorbate.bead_names[b])
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    
    # get epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)
    
    if adsorbate.nbeads == 1  # no need to perform rotations
        function E(fractional_coord::Array{Float64})
            """
            Energy (Boltzmann weighted at pos x
            """
            adsorbate.translate_to(framework.f_to_cartesian_mtrx * fractional_coord)
            return _vdW_energy_of_adsorbate!(adsorbate, epsilons, sigmas, framework, rep_factors, cutoff) 
        end
    else  # rotations necessary!
        function E(fractional_coord::Array{Float64})
            """
            Energy (Boltzmann weighted at pos x
            """
            adsorbate.translate_to(framework.f_to_cartesian_mtrx * fractional_coord)
            boltzmann_weight_sum = 0.0
            weighted_energy_sum = 0.0
            for k = 1:num_rotation_samples
                adsorbate.perform_uniform_random_rotation()
                _energy = _vdW_energy_of_adsorbate!(adsorbate, epsilons, sigmas, framework, rep_factors, cutoff) 
                boltzmann_weight = exp(-_energy / temperature)
                boltzmann_weight_sum += boltzmann_weight
                weighted_energy_sum += boltzmann_weight * _energy
            end
            return weighted_energy_sum / boltzmann_weight_sum
        end
    end
    
    res = optimize(E, x_f_start, method=:nelder_mead)
    
    # If optimization routine found a point outside of the home unit cell, try a new starting point
    while ((sum(res.minimum .< [-0.000001, -0.00001, -0.000001]) != 0) | (sum(res.minimum .> [1.00001, 1.000001, 1.000010]) != 0))
        for i = 1:3
            x_f_start[i] = mod(res.minimum[i], 1.0)
            x_f_start[i] += 0.05 * rand()
        end
        
        @printf("Fractional coords went outside of unit box; trying another starting point (%f, %f, %f)\n", x_f_start[1], x_f_start[2], x_f_start[3])
        res = optimize(E, x_f_start, method=:nelder_mead)
 #         error("Fractional coords went outside of unit box; choose better starting point\n")
        # bring into home unit cell (Fractional coords in [0,1]) and perturb randomly
    end
    @printf("Minimum E= %f at (%f, %f, %f)\n", res.f_minimum * 8.314/1000.0, res.minimum[1], res.minimum[2], res.minimum[3])

    return res.f_minimum * 8.314 / 1000.0, res.minimum  # minE, x_min
end

function get_optimal_rotation(adsorbate::Adsorbate, 
                              structurename::String, 
                              forcefieldname::String; 
                              num_rotation_samples::Int=500, 
                              cutoff::Float64=12.5)
    """
    Explore different rotations of the adsorbate and record the configuration of the one with the lowest energy
    """
    if (adsorbate.nbeads == 1)
        error("Dude, don't need rotations for a spherical adsorbate\n")
    end

    framework = Framework(structurename)
    
    @printf("Constructing forcefield object for %s...\n", forcefieldname)
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        @printf("\tBead %s...\n", adsorbate.bead_names[b])
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    
    # get epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)
    
    E_min = _vdW_energy_of_adsorbate!(adsorbate, epsilons, sigmas, framework, rep_factors, cutoff)
    opt_bead_xyz = adsorbate.bead_xyz
    for i = 1:num_rotation_samples
        adsorbate.perform_uniform_random_rotation()
        E = _vdW_energy_of_adsorbate!(adsorbate, epsilons, sigmas, framework, rep_factors, cutoff) 
        if E < E_min
            E_min = E
            opt_bead_xyz = adsorbate.bead_xyz
        end
    end

    adsorbate.bead_xyz = opt_bead_xyz
    
    return adsorbate, E_min
end
