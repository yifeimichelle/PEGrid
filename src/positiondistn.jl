###
#   Author: Cory M. Simon (corymsimon@gmail.com)
#   modified by Michelle Liu (meeechelle@gmail.com)
###
include("framework.jl")


function writeprobabilitydistncube(framework::Framework,
                                   adsorbatepositionfilename::AbstractString,
                                   which_adsorbate::AbstractString;
                                   binspacing::Float64=1.0,
                                   outputcubename::AbstractString="",
                                   molecule::AbstractString="",
                                   rep_factor::Int=1,
                                   z_symm::Bool=false,
                                   shift_x::Float64=0.0,
                                   shift_y::Float64=0.0,
                                   shift_z::Float64=0.0)
    """
    Takes an xyz file of adsorbate positions in the framework,
    bins these positions in 3D space
    outputs a Gaussian cube file of the probability of finding a particle in that voxel
    For make volume plots to visualize where molecules adsorb in the material
    """
    ### Partition unit cell into voxels
    # how many voxels in each direction?
    N_x = floor(Int, rep_factor * framework.a / binspacing) + 1
    N_y = floor(Int, rep_factor * framework.b / binspacing) + 1
    N_z = floor(Int, rep_factor * framework.c / binspacing) + 1
    @printf("Unit cell is partitioned into %d by %d by %d voxels, a total of %d voxels.\n", N_x, N_y, N_z, N_x*N_y*N_z)

    # spacing in fractional coords
    dx_f = rep_factor / N_x
    dy_f = rep_factor / N_y
    dz_f = rep_factor / N_z

    ### Store particle count in each voxel here
    counts = zeros(Float64, N_x, N_y, N_z)

    ### Loop through positions in .xyz file of adsorbate positions, store bin counts
    @printf("Reading xyz file of adsorbate positions, binning positions\n")
    if ~ isfile(adsorbatepositionfilename)
        @printf("Adsorbate positions .xyz file %s not present in working directory", adsorbatepositionfilename)
    end

    xyzfile = open(adsorbatepositionfilename, "r")
    readline(xyzfile)
    readline(xyzfile) # throw away two lines

    N_positions = 0
    for line in eachline(xyzfile)
        if length(split(line)) != 4
            continue
        end
        adsorbate = split(line)[1]

        if adsorbate == which_adsorbate
            # get (x, y, z) of adsorbate
            N_positions += 1
            x = float(split(line)[2]) - shift_y
            y = float(split(line)[3]) - shift_y
            z = float(split(line)[4]) - shift_z
            if z_symm
                z += framework.c*0.5
            end

            # get fractional coords
            f_coords = framework.cartesian_to_f_mtrx * [x, y, z]

            # reflect back to [0, rep_factor] for unit cell
            f_coords[1] = mod(f_coords[1], rep_factor)
            f_coords[2] = mod(f_coords[2], rep_factor)
            f_coords[3] = mod(f_coords[3], rep_factor)

            # get voxel index to which particle belongs
            i = floor(Int, 1 + f_coords[1] / dx_f)
            j = floor(Int, 1 + f_coords[2] / dy_f)
            k = floor(Int, 1 + f_coords[3] / dz_f)

            # check for particles out of bound due to numerical rounding
            if (i < 1)
                f_coords[1] += 0.001
            end
            if (j < 1)
                f_coords[2] += 0.001
            end
            if (k < 1)
                f_coords[3] += 0.001
            end
            if (i > N_x)
                f_coords[1] -= 0.001
            end
            if (j > N_y )
                f_coords[2] -= 0.001
            end
            if (k > N_z)
                f_coords[3] -= 0.001
            end

            # get voxel index to which particle belongs
            i = floor(Int, 1 + f_coords[1] / dx_f)
            j = floor(Int, 1 + f_coords[2] / dy_f)
            k = floor(Int, 1 + f_coords[3] / dz_f)

            if (i < 1)  |  (j < 1) | (k < 1)
                @printf("%s",line)
                @printf("%d %d %d\n", i, j, k)
                error("adsorbate outside of unit cell")
            end
            if (i > N_x)  |  (j > N_y ) | (k > N_z)
                @printf("%s",line)
                @printf("%d %d %d\n", i, j, k)
                error("adsorbate outside of unit cell")
            end

            # update counts
            counts[i, j, k] += 1
        end
    end

    @assert(sum(counts) == N_positions)

    println("Max number of adsorbates seen in a bin: ", maximum(counts))
    println("Fraction of bins that are empty: ", sum(counts .== 0) / length(counts))
    println("Number of non-empty bins: ", countnz(counts))
    println("Number of bins that have one adsorbate: ", sum(counts .== 1))
    # normalize to make a probability dist'n
    counts = counts / sum(counts)

    @printf("%d adsorbate %s positions recorded\n", N_positions, which_adsorbate);

    # Write cube file of probability distn
    if outputcubename == ""
        outputcubename = framework.structurename * "_" * molecule * which_adsorbate * "_adsorbate_probability_distn.cube"
    end
    probdistn = open(homedir() * "/PEGrid_output/" * outputcubename, "w")

    # Format of .cube described here http://paulbourke.net/dataformats/cube/
    write(probdistn, "This is an adsorbate probability distn file generated by PEviz\nLoop order: x, y, z\n")
    @printf(probdistn, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)  # 0 atoms, then origin
    # TODO list atoms in the crystal structure
    @printf(probdistn, "%d %f %f %f\n" , N_x, rep_factor * framework.f_to_cartesian_mtrx[1,1] / (N_x-1), 0.0, 0.0)  # N_x, vector along x-edge of voxel
    @printf(probdistn, "%d %f %f %f\n" , N_y, rep_factor * framework.f_to_cartesian_mtrx[1,2] / (N_y-1), rep_factor * framework.f_to_cartesian_mtrx[2,2] / (N_y-1), 0.0)  # N_y, vector along y-edge of voxel
    @printf(probdistn, "%d %f %f %f\n" , N_z, rep_factor * framework.f_to_cartesian_mtrx[1,3] / (N_z-1), rep_factor * framework.f_to_cartesian_mtrx[2,3] / (N_z-1), rep_factor * framework.f_to_cartesian_mtrx[3,3] / (N_z-1))

    @printf("Writing adsorbate probability distn to cube file:\n\t%s\n\t... ... ...\n", homedir() * "/PEGrid_output/" * outputcubename)
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
    @printf("\tDone. If visualizing in VisIt, do not uncheck the Extend unit cell box when opening the cube.\n")

end

function convert_pdb_to_xyz(inputfilename::AbstractString, outputfilename::AbstractString)
    """
    Takes a PDB file and converts it to XYZ format
    """
    inputfile = open(inputfilename)
    outputfile = open(outputfilename, "w")
    num_lines = readall(pipeline(`grep "ATOM" $inputfilename`, `wc -l`))
    @printf(outputfile, "%s\n", num_lines)
    for line in eachline(inputfile)
        s = split(line)
        if s[1] == "ATOM"
            @printf(outputfile, "%s %s %s %s\n", s[3], s[5], s[6], s[7])
        end
    end
    close(inputfile)
    close(outputfile)
    @printf("\tDone. Output written to %s\n", outputfilename)
end
