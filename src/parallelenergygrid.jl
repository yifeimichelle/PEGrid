function evaluate_expressions_on_all_cores(;kwargs...)
    """
    Evaluates a given expression on each process. For broadcasting variables
    """
    for (nm, val) in kwargs # unpacks expressions in args
        for p = 1:nprocs() #in workers()
            @spawnat p eval(Main, Expr(:(=), nm, val))
        end
    end
end

function parallel_writegrid(adsorbatename::AbstractString, 
                            structurename::AbstractString, 
                            forcefieldname::AbstractString; 
                            gridspacing=0.1, 
                            cutoff=12.5, 
                            temperature::Float64=-1.0, 
                            num_rotation_samples::Int=1000)
    """
    Compute the potential energy of an adsorbate molecule on a 3D grid of points superimposed on the unit cell of the structure.
    Parallelized across cores.

    The grid is written to a file `structurename.cube`, in Gaussian cube format. The units of the energy are kJ/mol.

    :param: AbstractString adsorbate: the name of the adsorbate molecule, corresponding to the forcefield file
    """
    @printf("Number of parallel cores: %d\n", nprocs())
    require("src/framework.jl")  # these statements load these files on all cores
    require("src/forcefield.jl")
    require("src/energyutils.jl")
    require("src/adsorbates.jl")

    # broadcast function arguments to all cores
    evaluate_expressions_on_all_cores(adsorbatename=adsorbatename)
    evaluate_expressions_on_all_cores(structurename=structurename)
    evaluate_expressions_on_all_cores(forcefieldname=forcefieldname)
    evaluate_expressions_on_all_cores(gridspacing=gridspacing)
    evaluate_expressions_on_all_cores(cutoff=cutoff)
    evaluate_expressions_on_all_cores(temperature=temperature)
    evaluate_expressions_on_all_cores(num_rotation_samples=num_rotation_samples)
  
    # load the adsorbate on all cores
    @printf("Constructing adsorbate %s on all cores...\n", adsorbatename)
    @everywhere begin
        adsorbate = Adsorbate(adsorbatename)
        if (adsorbate.nbeads > 1) & (temperature == -1.0)
            error("Provide temperature for Boltzmann weighted rotations, nbeads > 1 in adsorbate.")
        end
    end

    # load the framework on all cores
    @printf("Constructing framework object on all cores for %s...\n", structurename)
    @everywhere framework = Framework(structurename)

    # load the forcefield on all cores
    @printf("Constructing forcefield(s) for bead(s) in %s...\n", forcefieldname)
    @everywhere begin
        forcefields = Forcefield[]  # list of forcefields
        for b = 1:adsorbate.nbeads
            push!(forcefields, Forcefield(forcefieldname, adsorbate.bead_names[b], cutoff=cutoff))
        end
    end

    # get unit cell replication factors for periodic BCs on all cores
    @everywhere rep_factors = get_replication_factors(framework.f_to_cartesian_mtrx, cutoff)
    @printf("Unit cell replication factors for LJ cutoff of %.2f A: %d by %d by %d\n", cutoff, rep_factors[1], rep_factors[2], rep_factors[3])

    # get array of framework atom positions and corresponding epsilons and sigmas for speed on all cores
    # TODO is a shared array better? is this too much memory?
    @everywhere epsilons, sigmas = _generate_epsilons_sigmas(framework, forcefields)

    ### ------------------------------------------------------------------------- ###
    ### This is all done on current core only. e.g. writing to grid file

    # how many grid points in each direction? 
    N_x = int(framework.a / gridspacing) + 1
    @everywhere N_y = int(framework.b / gridspacing) + 1
    @everywhere N_z = int(framework.c / gridspacing) + 1
    @printf("Grid is %d by %d by %d points, a total of %d grid points.\n", N_x, N_y, N_z, N_x*N_y*N_z)

    # fractional grid point spacing. Think of grid points as center of voxels.
    dx_f = 1.0 / (N_x - 1)
    dy_f = 1.0 / (N_y - 1)
    dz_f = 1.0 / (N_z - 1)
    @printf("Fractional grid spacing: dx_f = %f, dy_f = %f, dz_f = %f\n", dx_f, dy_f, dz_f)

    # get fractional coords of energy grid. 
    xf_grid = linspace(0.0, 1.0, N_x)
    @everywhere {
        yf_grid = linspace(0.0, 1.0, N_y)
        zf_grid = linspace(0.0, 1.0, N_z)
    }

    # get grid point spacing in Cartesian space, just for kicks ^.^
    cartesian_spacing = framework.f_to_cartesian_mtrx * [xf_grid[2]-xf_grid[1], yf_grid[2]-yf_grid[1], zf_grid[2]-zf_grid[1]]
    @printf("Grid spacing: dx = %.2f, dy = %.2f, dz = %.2f\n", cartesian_spacing[1], cartesian_spacing[2], cartesian_spacing[3])

    # open grid file
    if ! isdir(homedir() * "/PEGrid_output/" * forcefieldname)
       mkdir(homedir() * "/PEGrid_output/" * forcefieldname) 
    end
    gridfilename = homedir() * "/PEGrid_output/" * forcefieldname * "/" * framework.structurename * "_" * adsorbatename * ".cube"
    gridfile = open(gridfilename, "w")

    # Format of .cube described here http://paulbourke.net/dataformats/cube/
    write(gridfile, "This is a grid file generated by PEviz\nLoop order: x, y, z\n")
    @printf(gridfile, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)  # 0 atoms, then origin
    # TODO list atoms in the crystal structure
    @printf(gridfile, "%d %f %f %f\n" , N_x, framework.f_to_cartesian_mtrx[1,1] / (N_x - 1), 0.0, 0.0)  # N_x, vector along x-edge of voxel
    @printf(gridfile, "%d %f %f %f\n" , N_y, framework.f_to_cartesian_mtrx[1,2] / (N_y - 1), framework.f_to_cartesian_mtrx[2,2] / (N_y - 1), 0.0)  # N_y, vector along y-edge of voxel
    @printf(gridfile, "%d %f %f %f\n" , N_z, framework.f_to_cartesian_mtrx[1,3] / (N_z - 1), framework.f_to_cartesian_mtrx[2,3] / (N_z - 1), framework.f_to_cartesian_mtrx[3,3] / (N_z - 1))

    ### ------------------------------------------------------------------------- ###

    @printf("Writing grid...\n")
    # loop over [fractional] grid points, compute energies
        
    # pre-allocate 2D array that stores sheet of energies in y-z plane
    @everywhere E_yz_sheet = zeros(N_y, N_z)

    @everywhere function compute_E_on_sheet(xf::Float64)
        """
        Compute energy on sheet of points.

        N_y, N_z, pos_array, epsilons, sigmas, E_yz_sheet assumed present
        """
        for j = 1:N_y  # loop over y_f-grid points
            for k = 1:N_z  # loop over z_f-grid points
                
                # translate adsorbate to grid pt
                adsorbate.translate_to(framework.f_to_cartesian_mtrx * [xf, yf_grid[j], zf_grid[k]])

                # compute energy here
                if adsorbate.nbeads == 1  # no need for sampling rotations
                    E_yz_sheet[j, k] = _vdW_energy_of_adsorbate!(adsorbate,
                                              epsilons,
                                              sigmas,
                                              framework,
                                              rep_factors,
                                              cutoff)
                else  # need to sample rotations and get boltzmann-weighted average
                    boltzmann_weight_sum = 0.0
                    weighted_energy_sum = 0.0
                    for rot = 1:num_rotation_samples
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
                    E_yz_sheet[j, k] = weighted_energy_sum / boltzmann_weight_sum
                end
            end
        end
        
        return E_yz_sheet
    end
        
    np = nprocs()  # determine the number of processes available

    i = 0 # count of x_f grid points completed
    while i < N_x  # loop over x_f-grid points
        @printf("\tPercent finished: %.1f\n", 100.0*i/N_x)  # print progress
       
        # map np sheets on the y-z plane to np cores
        sheets = pmap(compute_E_on_sheet, xf_grid[i+1:minimum([i + np, N_x])])
        
        i += np # update count of x_f grid points completed

        # write results to grid file
        for s = 1:size(sheets,1) # loop over sheets
            for j in 1:N_y  # loop over y_f-grid points
                for k in 1:N_z  # loop over z_f-grid points
                    # write energy at this point to grid file
                    @printf(gridfile, "%e ", sheets[s][j,k] * 8.314 / 1000.0)  # store in kJ/mol
                    if (k % 6) == 0
                        @printf(gridfile, "\n")
                    end
                end # end loop in z_f-grid points
                @printf(gridfile, "\n")  # new line after z loop
            end
        end

    end

    close(gridfile)
    @printf("\tDone.\nGrid available in %s\n", gridfilename)
end
