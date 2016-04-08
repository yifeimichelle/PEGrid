###
#  Author: Cory M. Simon (CoryMSimon@gmail.com)
###
include("framework.jl")
include("forcefield.jl")
include("adsorbates.jl")
include("getEwaldparams.jl")
using Optim
using ParallelAccelerator
    
# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const epsilon_0 = 4.7622424954949676e-7;  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)

function _electrostatic_potential(xyz_coord::Array{Float64},
                                 framework::Framework,
                                 rep_factors::Array{Int, 1},
                                 sr_cutoff::Float64,
                                 alpha::Float64,
                                 k_reps::Array{Int};
                                 verboseflag::Bool=false)
    """
    Compute electrostatic energy at a point inside a Framework using EWald summations.
    
    Parameters:
        xyz_coord: the Cartesian (x,y,z) coordinate inside framework where we seek electrostatic potential
        framework: the crystal structure in PEGrid Framework type
        rep_factors: x,y,z replication factors of primitive unit cell for PBCs for short-range energy calculations
        srcutoff: short-range cutoff distance
        alpha: decay parameter for long-range calcs
        k_reps: number of replications in k-space for long-range calcs
    """
    # get fractional coord of bead
    fractional_coord::Array{Float64, 1} = framework.cartesian_to_f_mtrx * xyz_coord
    
    # reflect bead to [0,1] via periodic boundary conditions
    fractional_coord[1] = mod(fractional_coord[1], 1.0)
    fractional_coord[2] = mod(fractional_coord[2], 1.0)
    fractional_coord[3] = mod(fractional_coord[3], 1.0)

    ######
    ######  Short-range
    ######
    const sqrt_alpha = sqrt(alpha) # compute sqrt(\alpha) here to save time
    E_sr::Float64 = 0.0  # short range energy
    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3]
                # get vector between this point and framework atoms in fractional space
                @inbounds dx::Array{Float64} = broadcast(-, fractional_coord + 1.0 * [rep_x, rep_y, rep_z], framework.fractional_coords)

                # transfer to Cartesian coords
                @inbounds dx = framework.f_to_cartesian_mtrx * dx

                @inbounds for i = 1:framework.natoms
                    @inbounds r::Float64 = norm(dx[:, i])
                    if r < sr_cutoff
                        @inbounds E_sr += framework.charges[i] / r * erfc(r * sqrt_alpha)
                    end
                end
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction
    E_sr /= 4.0 * pi * epsilon_0

    ######
    ######  Long-range
    ######
    # Cartesian vector from point to framework atoms (in matrix)
    @inbounds x_framework_to_pt::Array{Float64} = broadcast(-, xyz_coord, framework.f_to_cartesian_mtrx * framework.fractional_coords)
    E_lr::Float64 = 0.0
    # use symmetry: cos(x) = cos(-x). So go from 0:k_reps in x-direction, multiply by two except for kx == 0.
    # see here http://physics.stackexchange.com/questions/205794/symmetry-in-program-for-ewald-summation
    for kx = 0:k_reps[1]
        for ky = -k_reps[2]:k_reps[2]
            for kz = -k_reps[3]:k_reps[3]
                if (kx == 0) & (ky == 0) & (kz == 0)
                    # do not include home box
                    continue
                end
                
                # get the k-vector
                k::Array{Float64, 1} = framework.reciprocal_lattice * [kx, ky, kz]
                magnitude_k_2::Float64 = dot(k, k)  # its magnitude, squared
                
                # dot product of k-vector and dx
                @inbounds k_dot_dx::Array{Float64} = transpose(k) * x_framework_to_pt
                
                E_lr_this_k::Float64 = 0.0
                
                @simd for i = 1:framework.natoms
                    @inbounds E_lr_this_k += framework.charges[i] * cos(k_dot_dx[i])
                end
                
                # two is for symmetry!
                if (kx == 0)
                    E_lr_this_k = E_lr_this_k / magnitude_k_2 * exp(- magnitude_k_2 / (4.0 * alpha))
                else
                    E_lr_this_k = 2.0 * E_lr_this_k / magnitude_k_2 * exp(- magnitude_k_2 / (4.0 * alpha))
                end
                E_lr += E_lr_this_k
            end
        end
    end

    E_lr = E_lr / framework.v_unitcell / epsilon_0

    if verboseflag
        @printf("E, long range = %f K\n", E_lr)
        @printf("E, short range = %f K\n", E_sr)
    end
    return E_lr + E_sr
end

function _vdW_energy_of_bead(xyz_coord::Array{Float64},
                             epsilons::Array{Float64},
                             sigmas2::Array{Float64},
                             framework::Framework,
                             rep_factors::Array{Int, 1},
                             cutoff2::Float64)
    """
    Compute Van der Waals interaction energy via the LJ potential of a bead at xyz_coord (Cartesian!)
   
    Parameters:
        xyz_coord: Cartesian coordinate of bead. shape: (3,)
        epsilons: Lennard-Jones epsilon parameters for bead with corresponding framework atoms. shape: (framework.natoms,)
        sigmas2: Lennard-Jones sigma parameters (squared!) for bead with corresponding framework atoms. shape: (framework.natoms,)
        framework: the structure (PEGrid object)
        rep_factors: x,y,z replication factors of primitive unit cell for PBCs
        cutoff2: Lennard-Jones cutoff distance, squared
    """
    # get fractional coord of bead
    fractional_coord::Array{Float64,1} = framework.cartesian_to_f_mtrx * xyz_coord

    # reflect bead to [0,1] via periodic boundary conditions
    fractional_coord[1] = mod(fractional_coord[1], 1.0)
    fractional_coord[2] = mod(fractional_coord[2], 1.0)
    fractional_coord[3] = mod(fractional_coord[3], 1.0)

    energy::Float64 = 0.0  # initialize energy of this bead as 0 to subsequently add pairwise contributions
    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3]
                # get vectors between this point and the framework atom (perhaps shifted to ghost unit cell)
                dx::Array{Float64} = broadcast(-, framework.fractional_coords, fractional_coord - 1.0 * [rep_x, rep_y, rep_z])
                
                # convert to Cartesian coords
                dx = framework.f_to_cartesian_mtrx * dx
                
                # this is the distance, squared
                r2::Array{Float64, 1} = vec(sum(dx .* dx, 1))

                @simd for i = 1:framework.natoms
                    if r2[i] < cutoff2
                        @inbounds sig_ovr_r6::Float64 = sigmas2[i] / r2[i]
                        sig_ovr_r6 = sig_ovr_r6 * sig_ovr_r6 * sig_ovr_r6
                        @inbounds energy += 4.0 * epsilons[i] * sig_ovr_r6 * (sig_ovr_r6 - 1.0)
                    end
                end
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction

    return energy
end

function _vdW_energy_of_adsorbate(adsorbate::Adsorbate,
            epsilons::Array{Float64}, 
            sigmas2::Array{Float64},
            framework::Framework,
            rep_factors::Array{Int, 1}, 
            cutoff::Float64)
    """
    Compute Van der Waals interaction energy via pairwise Lennard Jones potentials.
    
    :param Adsorbate adsorbate: the molecule (includes position) that we want to calc the energy of
    :param Array{Float64} epsilons: Lennard-Jones epsilon parameters corresponding to each framework atom in the array
    :param Array{Float64} sigmas2: Lennard-Jones sigma parameters (squared) corresponding to each framework atom in the array
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: Lennard-Jones cutoff distance
    :returns: Float64 E: energy of adsorbate at (x_f,y_f,z_f) in Kelvin (well, same units as epsilon)
    """
    energy = 0.0
    for b = 1:adsorbate.nbeads  # for each bead in adsorbate
        energy += _vdW_energy_of_bead(adsorbate.x[:, b], epsilons[:, b], sigmas2[:, b], framework, rep_factors, cutoff*cutoff)
    end  # end loop over beads in adsorbate

    return energy  # in Kelvin
end

function get_replication_factors(framework::Framework, cutoff::Float64)
    """
    Get replication factors for the unit cell such that only directly adjacent ghost unit cells 
        are required for consideration in applying periodic boundary conditions (PBCs).
    That is, a particle in the "home" unit cell will not be within the cutoff distance of 
        any particles two ghost cells away from the home cell.
    This is done by ensuring that a sphere of radius cutoff can fit inside the unit cell
    
    :param Framework framework: the crystal structure Framework object
    :param Float64 cutoff: Lennard-Jones cutoff distance, beyond which VdW interactions are approximated as zero.
    """
    # unit cell vectors
    a = framework.f_to_cartesian_mtrx[:, 1]
    b = framework.f_to_cartesian_mtrx[:, 2]
    c = framework.f_to_cartesian_mtrx[:, 3]

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
        a = a + framework.f_to_cartesian_mtrx[:, 1]
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # b replication
    while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff / 2.0
        rep_factors[2] += 1
        b = b + framework.f_to_cartesian_mtrx[:, 2]
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    # c replication
    while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff / 2.0
        rep_factors[3] += 1
        c = c + framework.f_to_cartesian_mtrx[:, 3]
        c0 = [a b c] * [0.5, 0.5, 0.5] # update center
    end
    
    return rep_factors
end

function _generate_epsilons_sigmas2(framework::Framework, forcefields::Array{Forcefield})
    """
    For speeding up energy computations, we use arrays of LJ parameters corresponding to adsorbate-framework atom i interactions.

    Arrays `epsilons` and `sigmas2` contain corresponding epsilon and sigma^2 LJ parameters for each framework atom.
    Parameters:
        framework: PEGrid Framework object; order will correpsond to the atoms in this object.
        forcefields: array of PEGrid Forcefield objects; one for each bead of the adsorbate
    """
    # one Forcefield for each bead in the adsorbate
    n_beads = length(forcefields)
    # store LJ parameters here
    epsilons = Array(Float64, (framework.natoms, n_beads))
    sigmas2 = Array(Float64, (framework.natoms, n_beads))

    for ff = 1:n_beads
        for i = 1:framework.natoms
            if ! haskey(forcefields[ff].epsilon, framework.atoms[i])
                error(@sprintf("Atom %s not present in force field.", framework.atoms[i]))
            end
            epsilons[i, ff] = forcefields[ff].epsilon[framework.atoms[i]]
            sigmas2[i, ff] = forcefields[ff].sigma[framework.atoms[i]] ^ 2
        end
    end

    return epsilons, sigmas2
end

function vdW_energy_of_adsorbate_config(adsorbate::Adsorbate, framework::Framework, forcefieldname::AbstractString; cutoff::Float64=12.5)
    """
    Get van der Waals energy of adsorbate inside a framework.
    Parameters:
        adsorbate: PEGrid Adsorbate object (contains Cartesian coords of adsorbate)
        framework: PEGrid Framework object
        forcefieldname: name of force field to use for Lennard-Jones parameters
    """
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework, cutoff)
    
    # get position array and epsilons/sigmas for easy computation
    epsilons, sigmas2 = _generate_epsilons_sigmas2(framework, forcefields)
    return _vdW_energy_of_adsorbate(adsorbate,
                                epsilons,
                                sigmas2,
                                framework,
                                rep_factors, 
                                cutoff)
end

function electrostatic_energy_of_adsorbate(adsorbate::Adsorbate, framework::Framework; 
                                           sr_cutoff::Float64=12.5, 
                                           rep_factors::Array{Int}=[-1,-1,-1], 
                                           ewald_params::Dict=Dict())
    """
    Get electrostatic energy of adsorbate inside a framework.
    Parameters:
        adsorbate: PEGrid Adsorbate object (contains Cartesian coords of adsorbate)
        framework: PEGrid Framework object
        sr_cutoff: short range cutoff
    """
    if ! adsorbate.charged_flag
        return 0.0
    end
    # get unit cell replication factors for periodic BCs
    if rep_factors == [-1, -1, -1]
        rep_factors = get_replication_factors(framework, sr_cutoff)
    end
    if ! (haskey(ewald_params, "alpha") & haskey(ewald_params, "k_reps"))
        ewald_params = getEwaldparams(framework, sr_cutoff, 1e-6, verboseflag=true)
    end

    energy = 0.0
    for i = 1:adsorbate.nbeads
        if adsorbate.charges[i] == 0.0
            continue
        end
        potential_in_framework = _electrostatic_potential(adsorbate.x[:, i],
                                         framework, rep_factors, sr_cutoff,
                                         ewald_params["alpha"], ewald_params["k_reps"])
        energy += adsorbate.charges[i] * potential_in_framework
   end
   return energy
end

function find_min_energy_position(framework::Framework,
                                  forcefieldname::AbstractString, 
                                  adsorbatename::AbstractString,
                                  x_f_start::Array{Float64};
                                  cutoff::Float64=12.5,
                                  verboseflag::Bool=false)
    """
    Find minimum energy position of adsorbate in framework

    Parameters:
        framework: PEGrid Framework object that contains crystal structure info
        forcefieldname: name of force field to use for vdW interactions
        adsorbatename: name of adsorbate whose min energy position we seek
        x_f_start: fractional coordinates at which we initilize the optimization routine (COM of adsorbate)
        cutoff: vdW LJ cutoff and for short-range Ewald summations
        verboseflag: print progress?
    """
    if verboseflag
        @printf("Constructing adsorbate %s...\n", adsorbatename)
    end
    adsorbate = Adsorbate(adsorbatename)
    
    if verboseflag
        @printf("Constructing forcefield object for %s...\n", forcefieldname)
    end
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        if verboseflag
            @printf("\tBead %s...\n", adsorbate.bead_names[b])
        end
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    const rep_factors = get_replication_factors(framework, cutoff)
    if verboseflag
        @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    end
    
    # get Ewald params in case adsorbate is charged
    const ewald_params = getEwaldparams(framework, cutoff, 1e-6, verboseflag=true)
    
    # get epsilons/sigmas for easy computation
    const epsilons, sigmas2 = _generate_epsilons_sigmas2(framework, forcefields)
    
    # create functions for optimization routine
    if adsorbate.nbeads == 1  # no need to perform rotations
        function energy(fractional_coord::Array{Float64})
            """
            This function computes the energy (vdW + electrostatic) of adsorboate at given fractional coords
            """
            adsorbate.translate_to(framework.f_to_cartesian_mtrx * fractional_coord)
            energy_vdW = _vdW_energy_of_adsorbate(adsorbate, epsilons, sigmas2, framework, rep_factors, cutoff) 
            energy_electrostatic = electrostatic_energy_of_adsorbate(adsorbate, framework; sr_cutoff=cutoff, rep_factors=rep_factors, ewald_params=ewald_params)
            return energy_vdW + energy_electrostatic
        end
    else  # rotations necessary!
        function energy(x::Array{Float64})
            """
            This function computes the energy (vdW + electrostatic) of adsorboate at given fractional coords and angles

            x: first three elements are fractional coords, remaining 3 are Euler angles
            """
            adsorbate.translate_and_rotate_to(framework.f_to_cartesian_mtrx * x[1:3], x[4], x[5], x[6])
            energy_vdW = _vdW_energy_of_adsorbate(adsorbate, epsilons, sigmas2, framework, rep_factors, cutoff)
            energy_electrostatic = electrostatic_energy_of_adsorbate(adsorbate, framework; sr_cutoff=cutoff, rep_factors=rep_factors, ewald_params=ewald_params)
            return energy_vdW + energy_electrostatic
        end
    end
    
    # give starting point for angles if adsorbate has more than 1 bead.
    if adsorbate.nbeads > 1
        x_f_start = vcat(x_f_start, [rand() * pi, rand() * pi, rand() * pi])
    end

    res = optimize(energy, x_f_start, method=:l_bfgs)#method=:nelder_mead)
    
    # If optimization routine found a point outside of the home unit cell, try a new starting point
    while ((sum(res.minimum[1:3] .< [-0.000001, -0.00001, -0.000001]) != 0) | (sum(res.minimum[1:3] .> [1.00001, 1.000001, 1.000010]) != 0))
        for i = 1:3
            x_f_start[i] = mod(res.minimum[i], 1.0)
            x_f_start[i] += 0.05 * rand()
            x_f_start[i] = mod(x_f_start[i], 1.0)
        end

        if adsorbate.nbeads > 1
            x_f_start[4] = pi * rand()
            x_f_start[5] = pi * rand()
            x_f_start[6] = pi * rand()
        end
        
        if verboseflag
            @printf("Fractional coords went outside of unit box; trying another starting point (%f, %f, %f)\n", x_f_start[1], x_f_start[2], x_f_start[3])
        else
            @printf("Found optimal coords outside home box, retrying...\n")
        end
        res = optimize(energy, x_f_start, method=:l_bfgs)
 #         error("Fractional coords went outside of unit box; choose better starting point\n")
        # bring into home unit cell (Fractional coords in [0,1]) and perturb randomly
    end
    x_min = framework.f_to_cartesian_mtrx * res.minimum[1:3]
    if verboseflag
        @printf("Minimum E= %f kJ/mol at xf = (%f, %f, %f); x = (%f, %f, %f)\n", 
            res.f_minimum * 8.314 / 1000.0, 
            res.minimum[1], res.minimum[2], res.minimum[3],
            x_min[1], x_min[2], x_min[3])
    end
   
    if adsorbate.nbeads > 1 
        adsorbate.translate_and_rotate_to(x_min, res.minimum[4], res.minimum[5], res.minimum[6])
    else
        adsorbate.to(x_min)
    end

    return adsorbate, res.f_minimum
end
