###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
using DataFrames


type Forcefield
    """
    Store attributes of a forcefield
    """
    name::String

    # adsorbate molecule
    adsorbate::String

    # number of interactions in the forcefield
    ninteractions::Int

    # Lennard-Jones parameters for adsorbate - X interactions
    atoms::Array{String}
    epsilon::Array{Float64}
    sigma::Array{Float64}

    # cutoff radius for Lennard-Jones potential
    cutoff::Float64

    # mixing rules to get adsorbate-solid interactions
    mixingrules::String
    
    # constructor
    function Forcefield(name::String, 
                        adsorbate::String; 
                        cutoff::Float64=12.5, 
                        mixingrules="Lorenz-Berthelot")
        """
        Importing force field file, fill atrributes
        """
        forcefield = new()
        forcefield.name = name
        forcefield.adsorbate = adsorbate
        forcefield.cutoff = cutoff
        forcefield.mixingrules = mixingrules

        # read in pure X-X iteractions data
        if ~ isfile("data/forcefields/" * name * ".csv")
            @printf("Could not find file data/forcefields/%s.csv", name)
        end
        df = readtable("data/forcefields/" * name * ".csv", allowcomments=true)

        # get adsorbate epsilon and sigma
        if ~ (adsorbate in df[:atom])
            error(@sprintf("Adsorbate %s not present in forcefield!", adsorbate))
        end
        adsorbate_eps = df[df[:atom] .== adsorbate, :][:epsilon][1]
        adsorbate_sig = df[df[:atom] .== adsorbate, :][:sigma][1]
        
        # initialize arrays in forcefield
        forcefield.ninteractions = size(df, 1)
        forcefield.atoms = Array(String, forcefield.ninteractions)
        forcefield.epsilon = zeros(Float64, forcefield.ninteractions)
        forcefield.sigma = zeros(Float64, forcefield.ninteractions)

        # compute and store adsorbate - X sigma/epsilon LJ params
        for i = 1:forcefield.ninteractions
            forcefield.atoms[i] = df[:atom][i]
            if mixingrules == "Lorenz-Berthelot"
                forcefield.epsilon[i] = sqrt(adsorbate_eps * df[:epsilon][i])
                forcefield.sigma[i] = (adsorbate_sig + df[:sigma][i]) / 2.0
            elseif mixingrules == "Surface"
                forcefield.epsilon[i] = df[:epsilon][i]
                forcefield.sigma[i] = df[:sigma][i]
                error("These mixing rules are not implemented")
            end
        end
        
        return forcefield
        end
end
