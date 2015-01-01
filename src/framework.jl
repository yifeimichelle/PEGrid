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

function constructframework(structurename::String)
    """
    Construct and fill in information on framework type from a .cssr crystal structure file

    :param String structurename: name of structure. Will try to import file structurename.cssr
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

function replicate_cssr_to_xyz(frameworkname::String; rep_factor::Int=1)
    """
    Converts a .cssr crystal structure file to .xyz
    Replicates unit cell into rep_factor by rep_factor by rep_factor supercell.
    The primitive unit cell will be in the middle.

    :param String frameworkname: name of crystal structure
    :param Int rep_factor: number of times to replicate unit cell
    """
    framework = constructframework(frameworkname)
    xyz_file = open(framework.structurename * ".xyz", "w")
    @printf(xyz_file, "%d\n\n", framework.natoms * rep_factor^3)

    for a = 1:framework.natoms
        for rep_x = -rep_factor:rep_factor
            for rep_y = -rep_factor:rep_factor
                for rep_z = -rep_factor:rep_factor
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
end
