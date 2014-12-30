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

function E_vdw_at_point(structurename::String, forcefieldname::String, adsorbate::String, 
                        x_f::Float64, y_f::Float64, z_f::Float64; cutoff::Float64=12.5)
    """
    Compute energy at fractional grid point (x_f, y_f, z_f) of adsorbate inside framework.
    This is for a single use (not efficient for multiple calls).
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

function _E_vdw_at_point!(x_f::Float64, y_f::Float64, z_f::Float64,
            pos_array::Array{Float64}, epsilons::Array{Float64}, sigmas::Array{Float64},
            framework::Framework,
            rep_factors::Array{Int}, cutoff::Float64)
    """
    Compute Van der Waals interaction energy via the Lennard Jones potential.
    Vectorized: this function takes an Array{Float64} of :math:`r^2` distances from each framework atom

    :param Float64 x_f,y_f,z_f: fractional coordinates at which we seek to compute energy
    :param Array{Float64} pos_array: 3 by natoms array of fractional coords of framework atoms in primitive unit cell
    :param Array{Float64} epsilons: Lennard-Jones epsilon parameters corresponding to each framework atom in the array
    :param Array{Float64} sigmas: Lennard-Jones sigma parameters corresponding to each framework atom in the array
    :param Array{Int} rep_factors: x,y,z replication factors of primitive unit cell for PBCs
    :param Float64 cutoff: Lennard-Jones cutoff distance
    :returns: Float64 E: energy of adsorbate at (x_f,y_f,z_f) in Kelvin (well, same units as epsilon)
    """
    E = 0.0  # initialize energy at this grid point
    # loop over adjacent unit cells to implement periodic boundary conditions
    for rep_x = -rep_factors[1]:rep_factors[1]
        for rep_y = -rep_factors[2]:rep_factors[2]
            for rep_z = -rep_factors[3]:rep_factors[3] 

                # vector of grid pt. Can effectively move grid point instead of adding rep_factors in fractional coords to framework atoms
                x_gridpt = [x_f + 1.0 * rep_x, 
                            y_f + 1.0 * rep_y,
                            z_f + 1.0 * rep_z]
                
                # what follows is vectorized over the framework atoms in the primitive unit cell
                # subtract from each framework atom position the grid point
                dx = broadcast(-, pos_array, x_gridpt)
                
                # transfer to Cartesian coords
                dx = framework.f_to_cartesian_mtrx * dx
                
                # compute distance squared between grid point and each framework atom, r2
                r2 = sum(dx .* dx, 1)

                # select which interactions are within the cutoff distance
                idx_within_cutoff = r2 .< cutoff * cutoff

                # compute VdW energy with Lennard-Jones potential.
                sig_ovr_r6 = sigmas[idx_within_cutoff] .* sigmas[idx_within_cutoff] ./ r2[idx_within_cutoff]  # (sigma / r )^2
                sig_ovr_r6 = sig_ovr_r6 .* sig_ovr_r6 .* sig_ovr_r6
                E += sum(4.0 * epsilons[idx_within_cutoff] .* sig_ovr_r6 .* (sig_ovr_r6 - 1.0))
            end  # end replication in x-direction
        end  # end replication in y-direction
    end  # end replication in z-direction

    return E  # in Kelvin
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

function writegrid(adsorbate::String, structurename::String, forcefieldname::String; gridspacing=0.1, gridfilename=None, cutoff=12.5)
    """
    Compute the potential energy of an adsorbate molecule on a 3D grid of points superimposed on the unit cell of the structure.

    The grid is written to a file `structurename.cube`, in Gaussian cube format. The units of the energy are kJ/mol.

    :param: String adsorbate: the name of the adsorbate molecule, corresponding to the forcefield file
    """
    @printf("Constructing framework object for %s...\n", structurename)
    framework = constructframework(structurename)

    @printf("Constructing forcefield object for %s...\n", forcefieldname)
    forcefield = constructforcefield(forcefieldname, adsorbate, cutoff=cutoff)
    
    # get unit cell replication factors for periodic BCs
    rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    @printf("Unit cell replication factors for LJ cutoff of %.2f A: %d by %d by %d\n", forcefield.cutoff, rep_factors[1], rep_factors[2], rep_factors[3])
	# how many grid points in each direction? 
	N_x = int(framework.a / gridspacing) + 1
	N_y = int(framework.b / gridspacing) + 1
	N_z = int(framework.c / gridspacing) + 1
    @printf("Grid is %d by %d by %d points, a total of %d grid points.\n", N_x, N_y, N_z, N_x*N_y*N_z)
    # fractional grid point spacing. Think of grid points as center of voxels.
    dx_f = 1.0 / (N_x - 1)
    dy_f = 1.0 / (N_y - 1)
    dz_f = 1.0 / (N_z - 1)
    @printf("Fractional grid spacing: dx_f = %f, dy_f = %f, dz_f = %f\n", dx_f, dy_f, dz_f)

	# get fractional coords of energy grid. 
	xf_grid = linspace(0.0, 1.0, N_x)
	yf_grid = linspace(0.0, 1.0, N_y)
	zf_grid = linspace(0.0, 1.0, N_z)
    
    # get grid point spacing in Cartesian space, just for kicks ^.^
    cartesian_spacing = framework.f_to_cartesian_mtrx * [xf_grid[2]-xf_grid[1], yf_grid[2]-yf_grid[1], zf_grid[2]-zf_grid[1]]
    @printf("Grid spacing: dx = %.2f, dy = %.2f, dz = %.2f\n", cartesian_spacing[1], cartesian_spacing[2], cartesian_spacing[3])

    # get array of framework atom positions and corresponding epsilons and sigmas for speed
    pos_array, epsilons, sigmas = _generate_pos_array_epsilons_sigmas(framework, forcefield)
	
    # open grid file
    gridfile = open(homedir() * "/" * framework.structurename * "_" * forcefield.adsorbate * ".cube", "w")

    # Format of .cube described here http://paulbourke.net/dataformats/cube/
    write(gridfile, "This is a grid file generated by PEviz\nLoop order: x, y, z\n")
    @printf(gridfile, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)  # 0 atoms, then origin
    # TODO list atoms in the crystal structure
    @printf(gridfile, "%d %f %f %f\n" , N_x, framework.f_to_cartesian_mtrx[1,1] / (N_x - 1), 0.0, 0.0)  # N_x, vector along x-edge of voxel
    @printf(gridfile, "%d %f %f %f\n" , N_y, framework.f_to_cartesian_mtrx[1,2] / (N_y - 1), framework.f_to_cartesian_mtrx[2,2] / (N_y - 1), 0.0)  # N_y, vector along y-edge of voxel
    @printf(gridfile, "%d %f %f %f\n" , N_z, framework.f_to_cartesian_mtrx[1,3] / (N_z - 1), framework.f_to_cartesian_mtrx[2,3] / (N_z - 1), framework.f_to_cartesian_mtrx[3,3] / (N_z - 1))

    @printf("Writing grid...\n")
	# loop over [fractional] grid points, compute energies
	for i in 1:N_x  # loop over x_f-grid points
        # print progress
        if i % (int(N_x/10.0)) == 0
            @printf("\tPercent finished: %.1f\n", 100.0*i/N_x)
        end

		for j in 1:N_y  # loop over y_f-grid points
			for k in 1:N_z  # loop over z_f-grid points
                E = _E_vdw_at_point!(xf_grid[i], yf_grid[j], zf_grid[k], 
                                    pos_array, epsilons, sigmas, 
                                    framework,
                                    rep_factors, cutoff)
                
                # write to grid file
                @printf(gridfile, "%e ", E * 8.314 / 1000.0)  # store in kJ/mol
                if (k % 6) == 0
                    @printf(gridfile, "\n")
                end
			end # end loop in z_f-grid points
            @printf(gridfile, "\n")
		end # end loop in y_f-grid points
	end # end loop in x_f-grid points
    close(gridfile)
    @printf("\tDone.")
end
