###
#   Author: CoryMSimon@gmail.com
#   Adsorbate type and methods
###
using DataFrames

type Adsorbate
    name::String  # corresponds to name in force field
    
    nbeads::Int  # number LJ spheres that the adsorbate consists of

    # nbeads by 3 array with bead_positions
    # convention: first bead is 0, 0, 0
    bead_positions::Array{Float64}
    bead_names::Array{String}  # corresponds to name in force field

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
                adsorbate.bead_positions = zeros(adsorbate.nbeads, 3)
            end
            if i == 3
                for b in 1:adsorbate.nbeads
                    push!(adsorbate.bead_names, split(line)[1+b])
                end
            end
            if i > 5
                print(split(line, ","))
                adsorbate.bead_positions[i - 5, 1] = parsefloat(split(line, ",")[1])
                adsorbate.bead_positions[i - 5, 2] = parsefloat(split(line, ",")[2])
                adsorbate.bead_positions[i - 5, 3] = parsefloat(replace(split(line, ",")[3], "\n", ""))
            end
        end

        close(f)

        if (adsorbate.bead_positions[1,1] != 0.0) | (adsorbate.bead_positions[1,2] != 0.0) | (adsorbate.bead_positions[1,3] != 0.0)
            error("First bead must be at 0,0,0 by convention")
        end

        return adsorbate
    end
end
