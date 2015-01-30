###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
using DataFrames


type Framework
    """
    Stores info of a crystal structure
    """
    structurename::String
    # unit cell sizes (Angstrom)
    a::Float64
    b::Float64
    c::Float64

    # unit cell angles
    alpha::Float64
    beta::Float64
    gamma::Float64

    # volume of the unit cell (Angstrom^3)
    v_unitcell::Float64

    # atom count (in primitive unit cell)
    natoms::Int
    
    # constructor
    Framework() = new()

    # fractional coordinates
    xf::Array{Float64}
    yf::Array{Float64}
    zf::Array{Float64}

    # atom identites
    atoms::Array{String}

    # transformation matrix from fractional to cartesian
    f_to_cartesian_mtrx::Array{Float64}
end

function constructframework(structurename::String; check_coords=true)
    """
    Construct and fill in information on Framework type from a .cssr crystal structure file

    :param String structurename: name of structure. Will try to import file structurename.cssr
    :param Bool check_coords: check that fractional coords are in [0.,1]
    """
    # construct the framework object
    framework = Framework()
    framework.structurename = structurename

    # open crystal structure file
    if ~ isfile("data/structures/" * structurename * ".cssr")
        @printf("Crystal structure file %s not present in data/structures/", structurename * ".cssr")
    end
    f = open("data/structures/" * structurename * ".cssr")

    # get unit cell sizes
    line = split(readline(f)) # first line is (a, b, c)
    framework.a = float(line[1])
    framework.b = float(line[2])
    framework.c = float(line[3])

    # get unit cell angles. Convert to radians
    line = split(readline(f))
    framework.alpha = float(line[1]) * pi / 180.0
    framework.beta = float(line[2]) * pi / 180.0
    framework.gamma = float(line[3]) * pi / 180.0

    # write transformation matrix from fractional to cartesian coords
    v_unit = sqrt(1.0 - cos(framework.alpha)^2 - cos(framework.beta)^2 - cos(framework.gamma)^2 +
            2 * cos(framework.alpha) * cos(framework.beta) * cos(framework.gamma)) # volume of unit parallelpiped
    framework.v_unitcell = v_unit * framework.a * framework.b * framework.c  # volume of unit cell
    framework.f_to_cartesian_mtrx = Array(Float64, (3,3))


    framework.f_to_cartesian_mtrx[1,1] = framework.a
    framework.f_to_cartesian_mtrx[1,2] = framework.b * cos(framework.gamma)
    framework.f_to_cartesian_mtrx[1,3] = framework.c * cos(framework.beta)
    framework.f_to_cartesian_mtrx[2,1] = 0.0
    framework.f_to_cartesian_mtrx[2,2] = framework.b * sin(framework.gamma)
    framework.f_to_cartesian_mtrx[2,3] = framework.c * (cos(framework.alpha) - cos(framework.beta) * cos(framework.gamma)) / sin(framework.gamma)
    framework.f_to_cartesian_mtrx[3,1] = 0.0
    framework.f_to_cartesian_mtrx[3,2] = 0.0
    framework.f_to_cartesian_mtrx[3,3] = framework.c * v_unit / sin(framework.gamma)
    
    # get atom count, initialize arrays holding coords
    framework.natoms = int(split(readline(f))[1])
    framework.atoms = Array(String, framework.natoms)
    framework.xf = zeros(Float64, framework.natoms)  # fractional coordinates
    framework.yf = zeros(Float64, framework.natoms)
    framework.zf = zeros(Float64, framework.natoms)

    # read in atoms and fractional coordinates
    readline(f) # waste a line
    for a = 1:framework.natoms
        line = split(readline(f))

        framework.atoms[a] = line[2]

        framework.xf[a] = float(line[3]) % 1.0 # wrap to [0,1]
        framework.yf[a] = float(line[4]) % 1.0
        framework.zf[a] = float(line[5]) % 1.0
        
        if check_coords # assert fractional coords in [0,1]
            @assert ((framework.xf[a] >= 0.0) & (framework.xf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
            @assert ((framework.yf[a] >= 0.0) & (framework.yf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
            @assert ((framework.zf[a] >= 0.0) & (framework.zf[a] <= 1.0)) "Fraction coords not in [0,1]!\n"
        end
    end
    
    close(f) # close file

    return framework
end

function crystaldensity(framework::Framework)
    """
    Computes crystal density of a framework (kg/m3)
    """
    if ! isfile("data/atomicmasses.csv")
        print("Atomic masses file data/atomicmasses.csv not present")
    end
    df = readtable("data/atomicmasses.csv")

    mass = 0.0 # count mass of all atoms here
    for a = 1:framework.natoms
        if ~ (framework.atoms[a] in df[:atom])
            error(@sprintf("Framework atom %s not present in data/atomicmasses.cv file", framework.atoms[a]))
        end
        mass += df[df[:atom] .== framework.atoms[a], :][:mass][1]
    end
    
    return mass / framework.v_unitcell * 1660.53892  # --> kg/m3
end

function replicate_cssr_to_xyz(frameworkname::String; rep_factors::Array{Int}=[1, 1, 1])
    """
    Converts a .cssr crystal structure file to .xyz
    Replicates unit cell into rep_factor by rep_factor by rep_factor supercell.
    The primitive unit cell will be in the middle.

    :param String frameworkname: name of crystal structure
    :param Array{Float64} rep_factors: number of times to replicate unit cell
    """
    framework = constructframework(frameworkname)
    xyz_file = open(framework.structurename * ".xyz", "w")
    @printf(xyz_file, "%d\n\n", framework.natoms * (2*rep_factors[1]+1)*(2*rep_factors[2]+1)*(2*rep_factors[3]+1))

    for a = 1:framework.natoms
        for rep_x = -rep_factors[1]:rep_factors[1]
            for rep_y = -rep_factors[2]:rep_factors[2]
                for rep_z = -rep_factors[3]:rep_factors[3]
                    x_f = [framework.xf[a] + 1.0*rep_x,
                           framework.yf[a] + 1.0*rep_y,
                           framework.zf[a] + 1.0*rep_z]
                    xyz = framework.f_to_cartesian_mtrx * x_f
                    @printf(xyz_file, "%s %f %f %f\n", framework.atoms[a], xyz[1], xyz[2], xyz[3])
                end
            end
        end
    end
    close(xyz_file)
    @printf("Wrote file %s\n", framework.structurename * ".xyz")
end

function write_unitcell_boundary_vtk(frameworkname::String)
    """
    Write unit cell boundary as a .vtk file for visualizing the unit cell.
    """
    framework = constructframework(frameworkname)
    vtk_file = open(framework.structurename * ".vtk", "w")
    
    # write first lines
    @printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\nASCII\nDATASET POLYDATA\nPOINTS 8 double\n") 

    # write points on boundary of unit cell
    for i=0:1
        for j=0:1
            for k=0:1
            x = [i, j, k] # fractional coordinate
            cellpoint = framework.f_to_cartesian_mtrx * x
            @printf(vtk_file, "%.3f %.3f %.3f\n", cellpoint[1], cellpoint[2], cellpoint[3])
            end
        end
    end

    # define connections
    @printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")
    close(vtk_file)
    @printf("Wrote file %s\n", framework.structurename * ".vtk")
end

function shift_unitcell_in_cssr(frameworkname::String; shift::Array{Float64}=[0.5, 0.5, 0.5])
    """
    Shift the unit cell.

    :param String frameworkname: name of crystal structure
    """
    @assert sum((shift .< 1.0) & (shift .>=0 )) == 3
    framework = constructframework(frameworkname)
    new_cssr_file = open("data/structures/" * framework.structurename * "_shifted.cssr", "w")
    @printf(new_cssr_file, "%f %f %f\n%f %f %f\n%d 0\nshifted by PEGrid\n", 
        framework.a, framework.b, framework.c,
        framework.alpha*180.0/pi, framework.beta*180.0/pi, framework.gamma*180.0/pi,
        framework.natoms)
    
    count = 1 
    for a = 1:framework.natoms
        # consider the unit cells to the "left" that are now part of this UC
        for rep_x = -1:0
            for rep_y = -1:0
                for rep_z = -1:0
                    # fractional coordinate
                    x_f = [framework.xf[a] + 1.0*rep_x + shift[1],
                           framework.yf[a] + 1.0*rep_y + shift[2],
                           framework.zf[a] + 1.0*rep_z + shift[3]]
                    if sum((x_f .< 1.0) & (x_f .>= 0.0)) == 3
                        @printf(new_cssr_file, "%d %s %f %f %f\n", 
                            count, framework.atoms[a], 
                            x_f[1], x_f[2], x_f[3])
                        count += 1
                    end
                end
            end
        end
    end
    close(new_cssr_file)
    @printf("Wrote file data/structures/%s\n", framework.structurename * "_shifted.cssr")
end

function put_cssr_coords_in_0_1(frameworkname::String)
    """
    In an incorrect .cssr, reflect coords to [0,1]
    """
    framework = constructframework(frameworkname, check_coords=false)
    new_cssr_file = open("data/structures/" * framework.structurename * "_corrected.cssr", "w")
    @printf(new_cssr_file, "%f %f %f\n%f %f %f\n%d 0\ncorrected by PEviz\n", 
        framework.a, framework.b, framework.c,
        framework.alpha*180.0/pi, framework.beta*180.0/pi, framework.gamma*180.0/pi,
        framework.natoms)
    
    count = 1 
    for a = 1:framework.natoms
        @printf(new_cssr_file, "%d %s %f %f %f\n", 
            count, framework.atoms[a], 
            mod(framework.xf[a], 1.0), mod(framework.yf[a], 1.0),mod(framework.zf[a], 1.0))
        count += 1
    end
    close(new_cssr_file)
    @printf("Wrote correct file with coords in [0,1] in data/structures/%s\n", framework.structurename * "_corrected.cssr")
end
