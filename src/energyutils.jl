###
#  Author: Cory M. Simon (CoryMSimon@gmail.com)
###
include("framework.jl")
include("forcefield.jl")


function get_replication_factors(f_to_cartesian_mtrx::Array{Float64}, cutoff::Float64)
    """
    Get replication factors for the unit cell such that only directly adjacent ghost unit cells are required for consideration in applying periodic boundary conditions (PBCs).
    That is, a particle in the "home" unit cell will not be within the cutoff distance of any particles two ghost cells away from the home cell.
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

function _generate_pos_array_epsilons_sigmas(framework::Framework, forcefield::Forcefield)
    """
    For speeding up energy computations, we use arrays of framework atom positions and their corresponding LJ parameters.

    Each column in pos_array is [xf, yf, zf] for a framework atom. In Julia, arrays stored column-wise; 3 by natoms is speedy.
    Arrays `epsilons` and `sigmas` contain corresponding epsilon and sigma LJ parameters for each framework atom.
    """
    # framework atom positions
    pos_array = Array(Float64,(3,framework.natoms)) # stores fractional coords of the primitive unit cell of the framework
    pos_array[1,:] = framework.xf
    pos_array[2,:] = framework.yf
    pos_array[3,:] = framework.zf
    
    # array of LJ parameters that correspond to framework atoms
    epsilons = Array(Float64, (1, framework.natoms)) # make same shape as r2...
    sigmas = Array(Float64, (1, framework.natoms)) 
    for a = 1:framework.natoms
        if ~ (framework.atoms[a] in forcefield.atoms)
            error(@sprintf("Atom %s not present in force field.", framework.atoms[a]))
        end
        epsilons[a] = forcefield.epsilon[forcefield.atoms .== framework.atoms[a]][1]
        sigmas[a] = forcefield.sigma[forcefield.atoms .== framework.atoms[a]][1]
    end

    return  pos_array, epsilons, sigmas
end

function E_vdw_at_point(structurename::String, forcefieldname::String, adsorbate::String, 
                        x_f::Float64, y_f::Float64, z_f::Float64; cutoff::Float64=12.5)
    """
    Compute energy at fractional grid point (x_f, y_f, z_f) of adsorbate inside framework.
    This is for a casual use, as it is not efficient for multiple calls.
    Use internal function _E_vdw_at_point!(...) for writing energy grids.
    """
    @printf("Constructing framework object for %s...\n", structurename)
    framework = constructframework(structurename)

    @printf("Constructing forcefield object for %s...\n", forcefieldname)
    forcefield = constructforcefield(forcefieldname, adsorbate, cutoff=cutoff)
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    
    # get position array and epsilons/sigmas for easy computation
    pos_array, epsilons, sigmas = _generate_pos_array_epsilons_sigmas(framework, forcefield)

    return _E_vdw_at_point(x_f, y_f, z_f,
            pos_array, epsilons, sigmas, 
            framework,
            rep_factors, cutoff)
end
