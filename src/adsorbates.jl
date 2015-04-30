###
#   Author: CoryMSimon@gmail.com
#   Adsorbate type and methods
###
using DataFrames

type Adsorbate
    name::String  # corresponds to name in force field
    
    nbeads::Int  # number LJ spheres that the adsorbate consists of

    # 3 by nbeads array with bead_xyz, Cartesian coords of bead
    # convention: first bead is 0, 0, 0
    bead_xyz::Array{Float64}  # CARTESIAN (so irrespective of framework)
    bead_names::Array{String}  # corresponds to name in force field

    translate::Function  # translate adsorbate by Cartesian vector x
    rotate::Function  # perform a rotation of the adsorbate

    function Adsorbate(name::String)
        """
        Constructor
        """
        adsorbate = new()
        
        # open adsorbate file
        if ~ isfile("data/adsorbates/" * name * ".adsorbate")
            @printf("Adsorbate file %s.adsorbate not present in data/adsorbates/", name)
        end
        f = open("data/adsorbates/" * name * ".adsorbate")
        
        for (i, line) in enumerate(eachline(f))
            if i == 1
                adsorbate.name = split(line)[2]
            end
            if i == 2
                adsorbate.nbeads = parseint(split(line)[2])
                adsorbate.bead_names = String[]
                adsorbate.bead_xyz = zeros(3, adsorbate.nbeads)
            end
            if i == 3
                for b in 1:adsorbate.nbeads
                    push!(adsorbate.bead_names, split(line)[1+b])
                end
            end
            if i > 5
                adsorbate.bead_xyz[1, i - 5] = parsefloat(split(line, ",")[1])
                adsorbate.bead_xyz[2, i - 5] = parsefloat(split(line, ",")[2])
                adsorbate.bead_xyz[3, i - 5] = parsefloat(replace(split(line, ",")[3], "\n", ""))
            end
        end

        close(f)

        if (adsorbate.bead_xyz[1,1] != 0.0) | (adsorbate.bead_xyz[2,1] != 0.0) | (adsorbate.bead_xyz[3,1] != 0.0)
            error("First bead must be at 0,0,0 by convention")
        end

        adsorbate.translate = function (x::Array{Float64})
            """
            Translate adsorbate by Cartesian vector x
            (function of self)
            """
            adsorbate.bead_xyz = broadcast(+, adsorbate.bead_xyz, x)
        end 

        adsorbate.rotate = function (alpha::Float64, beta::Float64, gamma::Float64)
            """
            Rotation about Euler angles
            """
        end

        return adsorbate
    end
end
