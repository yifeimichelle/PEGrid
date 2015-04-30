###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
using DataFrames


type Forcefield
    """
    Stores attributes of a force field
    """
    name::String

    # bead inadsorbate molecule
    bead::String

    # number of interactions in the forcefield
    ninteractions::Int

    # Lennard-Jones parameters for bead - X interactions
    atoms::Array{String}
    epsilon::Array{Float64}
    sigma::Array{Float64}

    # cutoff radius for Lennard-Jones potential
    cutoff::Float64

    # mixing rules to get bead-solid interactions
    mixingrules::String
    
    # constructor
    function Forcefield(name::String, 
                        bead::String; 
                        cutoff::Float64=12.5,
                        mixingrules="Lorenz-Berthelot")
        """
        Importing force field file, fill atrributes
        """
        forcefield = new()
        forcefield.name = name
        forcefield.bead = bead
        forcefield.cutoff = cutoff
        forcefield.mixingrules = mixingrules

        # read in pure X-X iteractions data
        if ~ isfile("data/forcefields/" * name * ".csv")
            @printf("Could not find file data/forcefields/%s.csv", name)
        end
        df = readtable("data/forcefields/" * name * ".csv", allowcomments=true)

        # get bead epsilon and sigma
        if ~ (bead in df[:atom])
            error(@sprintf("Bead %s not present in forcefield!", bead))
        end
        bead_eps = df[df[:atom] .== bead, :][:epsilon][1]
        bead_sig = df[df[:atom] .== bead, :][:sigma][1]
        
        # initialize arrays in forcefield
        forcefield.ninteractions = size(df, 1)
        forcefield.atoms = Array(String, forcefield.ninteractions)
        forcefield.epsilon = zeros(Float64, forcefield.ninteractions)
        forcefield.sigma = zeros(Float64, forcefield.ninteractions)

        # compute and store bead - X sigma/epsilon LJ params
        for i = 1:forcefield.ninteractions
            forcefield.atoms[i] = df[:atom][i]
            if mixingrules == "Lorenz-Berthelot"
                forcefield.epsilon[i] = sqrt(bead_eps * df[:epsilon][i])
                forcefield.sigma[i] = (bead_sig + df[:sigma][i]) / 2.0
            elseif mixingrules == "PureInteractions"
                forcefield.epsilon[i] = df[:epsilon][i]
                forcefield.sigma[i] = df[:sigma][i]
            else
                error("These mixing rules are not implemented")
            end
        end
        
        return forcefield
    end
end
