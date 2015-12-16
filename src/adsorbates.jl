###
#   Author: CoryMSimon@gmail.com
#   Adsorbate type and methods
###
using DataFrames

function uniform_unit_vector_on_sphere()
    """
    Generate a uniformly distributed point on a sphere
    """
    v = zeros(3)
    while norm(v) < 0.0001
        v = [randn(), randn(), randn()]
    end
    return v / norm(v)
end

type Adsorbate
    name::AbstractString  # corresponds to name in force field
    
    nbeads::Int  # number LJ spheres that the adsorbate consists of

    # 3 by nbeads array with bead_xyz, Cartesian coords of bead
    # convention: first bead is 0, 0, 0
    bead_xyz::Array{Float64}  # CARTESIAN (so irrespective of framework). shape = (3, nbeads)
    bead_xyz_COM_origin::Array{Float64}  # Cartesian coords where origin is center of mass
    bead_names::Array{AbstractString}  # corresponds to name in force field
    COM::Array{Float64}  # store center of mass so we don't need to compute it each time

    translate_to::Function  # translate adsorbate by Cartesian vector x
    perform_uniform_random_rotation::Function  # perform a rotation of the adsorbate
    write_to_xyz::Function  # write adsorbate positions to .xyz

    get_MW::Function
    _get_COM::Function
    translate_COM_to_origin::Function  # set origin of coords at center of mass

    function Adsorbate(name::AbstractString)
        """
        Constructor
        """
        adsorbate = new()
        
        # open adsorbate file
        if ~ isfile("data/adsorbates/" * name * ".adsorbate")
            error(@sprintf("Adsorbate file %s.adsorbate not present in data/adsorbates/", name))
        end
        f = open("data/adsorbates/" * name * ".adsorbate")
        
        for (i, line) in enumerate(eachline(f))
            if i == 1
                adsorbate.name = split(line)[2]
            end
            if i == 2
                adsorbate.nbeads = parse(Int, split(line)[2])
                adsorbate.bead_names = AbstractString[]
                adsorbate.bead_xyz = zeros(3, adsorbate.nbeads)
            end
            if i == 3
                for b in 1:adsorbate.nbeads
                    push!(adsorbate.bead_names, split(line)[1+b])
                end
            end
            if i > 5
                adsorbate.bead_xyz[1, i - 5] = parse(Float64, split(line, ",")[1])
                adsorbate.bead_xyz[2, i - 5] = parse(Float64, split(line, ",")[2])
                adsorbate.bead_xyz[3, i - 5] = parse(Float64, replace(split(line, ",")[3], "\n", ""))
            end
        end

        close(f)

        if (adsorbate.bead_xyz[1,1] != 0.0) | (adsorbate.bead_xyz[2,1] != 0.0) | (adsorbate.bead_xyz[3,1] != 0.0)
            error("First bead must be at 0,0,0 by convention")
        end

        adsorbate.translate_to = function (x::Array{Float64})
            """
            Translate adsorbate from origin (at COM) to Cartesian vector x
            (function of self)
            """
            adsorbate.bead_xyz = broadcast(+, adsorbate.bead_xyz_COM_origin, x)
            adsorbate.COM = x  # update center of mass
        end 

        adsorbate.get_MW = function()
            """
            Get MW of adsorbate
            """
            if ! isfile("data/atomicmasses.csv")
                error("Atomic masses file data/atomicmasses.csv not present")
            end
            df_masses = readtable("data/atomicmasses.csv")
            
            MW = 0.0
            for i = 1:adsorbate.nbeads
                if ~ (adsorbate.bead_names[i] in df_masses[:atom])
                    error(@sprintf("Adsorbate bead %s not present in data/atomicmasses.cv file", 
                                adsorbate.bead_names[i]))
                end
               
                MW += df_masses[df_masses[:atom] .== adsorbate.bead_names[i], :][:mass][1]
            end
            return MW
        end

        adsorbate._get_COM = function()
            """
            Get center of mass of adsorbate
            """
            if ! isfile("data/atomicmasses.csv")
                error("Atomic masses file data/atomicmasses.csv not present")
            end
            df_masses = readtable("data/atomicmasses.csv")
            
            M = 0.0  # total mass
            x_COM = [0.0, 0.0, 0.0]
            for i = 1:adsorbate.nbeads
                if ~ (adsorbate.bead_names[i] in df_masses[:atom])
                    error(@sprintf("Adsorbate bead %s not present in data/atomicmasses.cv file", 
                                adsorbate.bead_names[i]))
                end
                
                m_bead = df_masses[df_masses[:atom] .== adsorbate.bead_names[i], :][:mass][1]
                M += m_bead
                x_COM += m_bead * adsorbate.bead_xyz[:, i]
            end
            return x_COM / M
        end
        
        # store bead coords with COM at origin as an attribute by subtracting COM
        adsorbate.bead_xyz_COM_origin = broadcast(-, adsorbate.bead_xyz, adsorbate._get_COM())

        adsorbate.translate_COM_to_origin = function()
            """
            Set bead_xyz at the center of mass
            """
            adsorbate.bead_xyz = adsorbate.bead_xyz_COM_origin
 #             eps = .00001
 #             assert(adsorbate._get_COM()[1] < eps)
 #             assert(adsorbate._get_COM()[2] < eps)
 #             assert(adsorbate._get_COM()[3] < eps)
 #             assert(adsorbate._get_COM()[1] > -eps)
 #             assert(adsorbate._get_COM()[2] > -eps)
 #             assert(adsorbate._get_COM()[3] > -eps)
        end
        
        adsorbate.perform_uniform_random_rotation = function ()
            """
            Perform a uniform random rotation of adsorbate (updating bead_xyz)
            Keep current center of mass
            # see http://www.mech.utah.edu/~brannon/public/rotation.pdf pg 106
            """
            # build rotation matrix R
            R = zeros(3,3)
            R[:, 1] = uniform_unit_vector_on_sphere()
            m = uniform_unit_vector_on_sphere()
            R[:, 2] = m - dot(m, R[:, 1]) * R[:, 1]  # subtract of component along col1 so it is orthogonal
            R[:, 2] = R[:, 2] / norm(R[:, 2])
            R[:, 3] = cross(R[:, 1], R[:, 2])  # gives orthogonal vector to first two
            
            adsorbate.bead_xyz = broadcast(+,  R * adsorbate.bead_xyz_COM_origin, adsorbate.COM)
        end

        adsorbate.write_to_xyz = function (filename::AbstractString)
            """
            Write adsorbate to .xyz file
            """
            filename = filename * ".xyz"

            f = open(filename, "w")
            write(f, @sprintf("%d\n\n", adsorbate.nbeads))
            for i = 1:adsorbate.nbeads
                write(f, @sprintf("%s %f %f %f\n", adsorbate.bead_names[i], adsorbate.bead_xyz[1, i], adsorbate.bead_xyz[2, i], adsorbate.bead_xyz[3, i]))
            end
            close(f)
        end

        adsorbate.translate_COM_to_origin()
        adsorbate.COM = adsorbate._get_COM()

        return adsorbate
    end
end
