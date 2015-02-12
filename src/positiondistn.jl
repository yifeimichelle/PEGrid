###
#   Author: Cory M. Simon (corymsimon@gmail.com)
###
include("framework.jl")


function writeprobabilitydistncube(structurename::String, which_adsorbate::String, binspacing=1)
    """
    Takes an xyz file of adsorbate positions in the framework,
    bins these positions in 3D space
    outputs a Gaussian cube file of the probability of finding a particle in that voxel
    For make volume plots to visualize where molecules adsorb in the material
    """
    ### load framework info
    framework = constructframework(structurename)

    ### Partition unit cell into voxels
	# how many points in each direction? 
	N_x = int(framework.a / binspacing) + 1
	N_y = int(framework.b / binspacing) + 1
	N_z = int(framework.c / binspacing) + 1
    @printf("Position disnt is %d by %d by %d points, a total of %d points.\n", N_x, N_y, N_z, N_x*N_y*N_z)

    # spacing in fractional coords
    dx_f = 1.0 / (N_x - 1)
    dy_f = 1.0 / (N_y - 1)
    dz_f = 1.0 / (N_z - 1)

    ### Store particle count in each voxel here
    counts = Array(Float64, N_x, N_y, N_z)

    ### Loop through positions in .xyz file of adsorbate positions, store bin counts
    @printf("Reading xyz file of adsorbate positions, binning positions\n")
    if ~ isfile("adsorbate_positions_" * structurename * ".xyz")
        @printf("Adsorbate positions .xyz file %s not present", "adsorbate_positions_" * structurename * ".xyz")
    end

    xyzfile = open("adsorbate_positions_" * structurename * ".xyz", "r")
    readline(xyzfile)
    readline(xyzfile) # throw away two lines
    
    N_positions = 0
    for line in eachline(xyzfile)
        adsorbate = split(line)[1]
        
        if adsorbate == which_adsorbate
            # get (x, y, z) of adsorbate
            N_positions += 1
            x = float(split(line)[2])
            y = float(split(line)[3])
            z = float(split(line)[4])
           
            # get fractional coords
            f_coords = framework.cartesian_to_f_mtrx * [x, y, z]

            # reflect back to [0, 1] for unit cell
            f_coords[1] = mod(f_coords[1], 1)
            f_coords[2] = mod(f_coords[2], 1)
            f_coords[3] = mod(f_coords[3], 1)

            # get bin index to which particle belongs
            # for Periodic Boundary Conditions, bins on edges are split to opposite faces
            i = 1 + floor(1 + (f_coords[1] - dx_f/2) / dx_f)
            j = 1 + floor(1 + (f_coords[2] - dy_f/2) / dy_f)
            k = 1 + floor(1 + (f_coords[3] - dz_f/2) / dz_f)

            # only fill first bin. copy to last bin later
            if i == N_x
                i = 1
            end
            if j == N_y
                j = 1
            end
            if k == N_z
                k = 1
            end

            if (i < 1)  |  (j < 1) | (k < 1)
                error("adsorbate outside of unit cell")
            end
            if (i > N_x)  |  (j > N_y ) | (k > N_z)
                error("adsorbate outside of unit cell")
            end
            
            # update counts
            counts[i, j, k] += 1
            
        end
    end

    # copy values in first point to last point since these are the same bin via PBC
    counts[N_x, :, :] = counts[1, :, :]
    counts[:, N_y, :] = counts[:, 1, :]
    counts[:, :, N_z] = counts[:, :, 1]

    counts = counts / sum(counts) # normalize

    @printf("%d adsorbate %s positions recorded\n", N_positions, which_adsorbate);

    # Write cube file of probability distn
    probdistn = open(homedir() * "/PEGrid_output/" * framework.structurename * "_adsorbate_probability_distn" * ".cube", "w")

    # Format of .cube described here http://paulbourke.net/dataformats/cube/
    write(probdistn, "This is an adsorbate probability distn file generated by PEviz\nLoop order: x, y, z\n")
    @printf(probdistn, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)  # 0 atoms, then origin
    # TODO list atoms in the crystal structure
    @printf(probdistn, "%d %f %f %f\n" , N_x, framework.f_to_cartesian_mtrx[1,1] / (N_x - 1), 0.0, 0.0)  # N_x, vector along x-edge of voxel
    @printf(probdistn, "%d %f %f %f\n" , N_y, framework.f_to_cartesian_mtrx[1,2] / (N_y - 1), framework.f_to_cartesian_mtrx[2,2] / (N_y - 1), 0.0)  # N_y, vector along y-edge of voxel
    @printf(probdistn, "%d %f %f %f\n" , N_z, framework.f_to_cartesian_mtrx[1,3] / (N_z - 1), framework.f_to_cartesian_mtrx[2,3] / (N_z - 1), framework.f_to_cartesian_mtrx[3,3] / (N_z - 1))

    @printf("Writing adsorbate probability distn to cube file:\n\t%s\n\t... ... ...\n", homedir() * "/PEGrid_output/" * framework.structurename * "_adsorbate_probability_distn" * ".cube")
	for i in 1:N_x  # loop over x_f-grid points
		for j in 1:N_y  # loop over y_f-grid points
			for k in 1:N_z  # loop over z_f-grid points

                # write probability at this point to grid file
                @printf(probdistn, "%e ", counts[i, j, k])
                if (k % 6) == 0
                    @printf(probdistn, "\n")
                end

			end # end loop in z_f-grid points
            @printf(probdistn, "\n")
		end # end loop in y_f-grid points
	end # end loop in x_f-grid points
    close(probdistn)
    @printf("\tDone.\n")
    
end