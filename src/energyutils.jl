###
#  Author: Cory M. Simon (CoryMSimon@gmail.com)
###
include("framework.jl")
include("forcefield.jl")
include("adsorbates.jl")


function _energy_of_bead(xyz_coord::Array{Float64},
                         epsilons::Array{Float64},
                         sigmas::Array{Float64},
                         framework::Framework,
                         rep_factors::Array{Int},
                         cutoff::Float64)
    """
    Compute Van der Waals interaction energy via the LJ potential of a bead at fractional_coord

    :param Array{Float64} epsilons: Lennard-Jones epsilon parameters for bead with corresponding framework atoms
    :param Array{Float64} sigmas: Lennard-Jones sigma parameters for bead with corresponding framework atoms
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: Lennard-Jones cutoff distance
    :returns: Float64 E: energy of adsorbate at (x_f,y_f,z_f) in Kelvin (well, same units as epsilon)
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

function _energy_of_adsorbate!(adsorbate::Adsorbate,
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
        E += _energy_of_bead(adsorbate.bead_xyz[:, b], epsilons[:, b], sigmas[:, b], framework, rep_factors, cutoff)
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
            if ~ (framework.atoms[i] in forcefields[ff].atoms)
                error(@sprintf("Atom %s not present in force field.", framework.atoms[i]))
            end
            epsilons[i, ff] = forcefields[ff].epsilon[forcefields[ff].atoms .== framework.atoms[i]][1]
            sigmas[i, ff] = forcefields[ff].sigma[forcefields[ff].atoms .== framework.atoms[i]][1]
        end
    end

    return  epsilons, sigmas
end

function energy_of_adsorbate(adsorbatename::String,
                        fractional_translations::Array{Float64},
                        structurename::String,
                        forcefieldname::String;
                        cutoff::Float64=12.5)
    """
    Compute energy of adsorbate

    :returns Array{Float64} of potential energies at each point in Kelvin (well, same units as \epsilon given in forcefield)
    """
    @printf("Constructing adsorbate %s...\n", adsorbatename)
    adsorbate = Adsorbate(adsorbatename)

    @printf("Constructing framework object for %s...\n", structurename)
    framework = Framework(structurename)

    @printf("Constructing forcefield(s) for bead(s) in %s...\n", forcefieldname)
    forcefields = Forcefield[]  # list of forcefields
    for b = 1:adsorbate.nbeads
        @printf("\tBead %s...\n", adsorbate.bead_names[b])
        push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
    end
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
    
    # get position array and epsilons/sigmas for easy computation
    epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)

    if (size(fractional_translations)[1] != 3)
        error("fractional_translations mtrx shld be 3 by npoints")
    end
    npoints = (length(fractional_translations) != 3) ? size(fractional_translations)[2] : 1

    @printf("\tComputing potential energy at %d points\n", npoints)

    E = zeros(npoints)  # pre-allocate array of energies

    for i = 1:npoints
        # translate adsorbate by vector:
        xyz_translation = framework.f_to_cartesian_mtrx * fractional_translations[:, i]
        adsorbate.translate(xyz_translation)
        # compute energy
        E[i] = _energy_of_adsorbate!(adsorbate,
                                    epsilons,
                                    sigmas,
                                    framework,
                                    rep_factors, 
                                    cutoff)

        # translate back
        adsorbate.translate(-xyz_translation)
        # TODO: cld numerical issues arise, subtracting adding many times? drift of orig coords...
    end
    
    return E  # in (Kelvin)
end

 # function find_min_energy_position(structurename::String, 
 #                                   forcefieldname::String, 
 #                                   adsorbatename::String,
 #                                   x_f_start::Array{Float64};
 #                                   cutoff::Float64=12.5)
 #     """
 #     Find minimum energy position of adsorbate in framework
 #     """
 #     @printf("Constructing framework object for %s...\n", structurename)
 #     framework = Framework(structurename)
 # 
 #     @printf("Constructing forcefield object for %s...\n", forcefieldname)
 #     forcefield = Forcefield(forcefieldname, adsorbate, cutoff=cutoff)
 # 
 #     @printf("Constructing adsorbate %s...\n", adsorbatename)
 #     adsorbate = Adsorbate(adsorbatename)
 #     if adsorbate.nbeads > 1
 #         error("Sorry not implemented for >1 beads yet...")
 #     end
 #     
 #     # get unit cell replication factors for periodic BCs
 #     rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
 #     @printf("\tUnit cell replication factors for cutoff radius %f A: %d x %d x %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
 #     
 #     # get epsilons/sigmas for easy computation
 #     epsilons, sigmas = _generate_epsilons_sigmas(framework, [forcefield])
 # 
 #     # define energy as a function of fractional coord
 #     E(x) = _E_vdw_at_point!(x,
 #                 epsilons, sigmas, 
 #                 framework,
 #                 rep_factors, cutoff)
 #     
 #     res = optimize(E, x_f_start, method=:nelder_mead)
 #     @printf("Minimum E= %f at (%f, %f, %f)\n", res.f_minimum * 8.314/1000.0, res.minimum[1], res.minimum[2], res.minimum[3])
 #     
 #     # If optimization routine found a point outside of the home unit cell, try a new starting point
 #     while ((sum(res.minimum .< [-0.000001, -0.00001, -0.000001]) != 0) | (sum(res.minimum .> [1.00001, 1.000001, 1.000010]) != 0))
 #         for i = 1:3
 #             x_f_start[i] = mod(res.minimum[i], 1.0)
 #             x_f_start[i] += 0.1 * rand()
 #         end
 #         
 #         @printf("Fractional coords went outside of unit box; trying another starting point (%f, %f, %f)\n", x_f_start[1], x_f_start[2], x_f_start[3])
 #         res = optimize(E, x_f_start, method=:nelder_mead)
 #  #         error("Fractional coords went outside of unit box; choose better starting point\n")
 #         # bring into home unit cell (Fractional coords in [0,1]) and perturb randomly
 #     end
 # 
 #     return res.f_minimum * 8.314 / 1000.0, res.minimum  # minE, x_min
 # end
