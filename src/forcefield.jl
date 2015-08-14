###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
using DataFrames


type Forcefield
    """
    Stores attributes of a Lennard-Jones force field
    """
    # name of force field
    name::String

    # bead in adsorbate molecule
    bead::String

    # number of interactions in the forcefield
    ninteractions::Int

    # Lennard-Jones parameters for bead - X interactions
    # store as dictionary so we can call e.g.:
    # epsilon["C"]
    # to get bead-C epsilon parameter
    epsilon::Dict
    sigma::Dict

    # cutoff radius for Lennard-Jones potential
    cutoff::Float64

    # mixing rules to get bead-solid interactions
    mixingrules::String
    
    print_info::Function
    
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
        forcefield.epsilon = Dict()
        forcefield.sigma = Dict()

        # compute and store bead - X sigma/epsilon LJ params
        for i = 1:forcefield.ninteractions
            atom_type = df[:atom][i]
            if mixingrules == "Lorenz-Berthelot"
                forcefield.epsilon[atom_type] = sqrt(bead_eps * df[:epsilon][i])
                forcefield.sigma[atom_type] = (bead_sig + df[:sigma][i]) / 2.0
            elseif mixingrules == "PureInteractions"
                forcefield.epsilon[atom_type] = df[:epsilon][i]
                forcefield.sigma[atom_type] = df[:sigma][i]
            else
                error("These mixing rules are not implemented")
            end
        end

        forcefield.print_info = function()
            @printf("%s Force Field for %s adsorbate bead.\n", forcefield.name, forcefield.bead)
            @printf("\tMixing rules: %s\n", forcefield.mixingrules)
            @printf("\tLJ cutoff radius: %f\n", forcefield.cutoff)

            for atom_type in keys(forcefield.epsilon)
                @printf("%-4s - %-4s. epsilon = %f K, sigma = %f A\n",
                    forcefield.bead, atom_type,
                    forcefield.epsilon[atom_type], forcefield.sigma[atom_type])
            end
        end
        
        return forcefield
    end
end
