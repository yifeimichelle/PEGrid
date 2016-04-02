###
#   Author: CoryMSimon@gmail.com
#   Adsorbate type and methods
###
using DataFrames
using Base.Test

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

function rotation_matrix(theta::Float64, direction::Int)
    """
    Return a matrix for a rotation about the x-, y-, or z- axis

    Parameters:
        theta: angle at which we'd like to rotate. degrees
        direction: 1, 2, or 3 for x, y, z axis.
    """
    theta = theta * pi / 180.0  # convert to radians
    R = zeros(Float64, 3, 3)
    cosine_ = cos(theta)
    sine_ = sin(theta)

    if direction == 1
        R[1, 1] = 1.0
        R[2, 2] = cosine_
        R[2, 3] = -sine_
        R[3, 2] = sine_
        R[3, 3] = cosine_
        return R
    end

    if direction == 2
        R[1, 1] = cosine_
        R[1, 3] = sine_
        R[2, 2] = 1.0
        R[3, 1] = -sine_
        R[3, 3] = cosine_
        return R
    end

    if direction == 3
        R[1, 1] = cosine_
        R[1, 2] = -sine_
        R[2, 1] = sine_
        R[2, 2] = cosine_
        R[3, 3] = 1.0
        return R
    end

    # if made it this far, didn't pass 1, 2, or 3 for the direction!
    error("Pass 1, 2, or 3 for the direction for x-, y-, and z- direction.\n")
end

type Adsorbate
    name::AbstractString  # corresponds to name in force field
    
    nbeads::Int  # number of Lennard-Jones beads in adsorbate
    charged_flag::Bool  # does the adsorbate have point charges?

    # 3 by nbeads array with bead_xyz, Cartesian coords of bead
    # convention: first bead is 0, 0, 0
    bead_xyz::Array{Float64}  # CARTESIAN (so irrespective of framework). shape = (3, nbeads)
    bead_names::Array{AbstractString}  # corresponds to name in force field
    bead_charges::Array{Float64}
    x_center_of_mass::Array{Float64}  # store center of mass so we don't need to compute it each time

    translate_to::Function  # translate adsorbate such that its center of mass is at x
    translate_by::Function # translate adsorbate by dx
    perform_uniform_random_rotation::Function  # perform a rotation of the adsorbate
    write_to_xyz::Function  # write adsorbate positions to .xyz
    rotate_about_axis::Function  # rotate adsorbate about a given axis
    print_bond_lengths::Function

    get_MW::Function
    _get_COM::Function
    center_at_origin ::Function  # subtract center of mass so that it is at origin

    function Adsorbate(name::AbstractString)
        """
        Constructor. Information for adsorbate is stored in file:
            data/adsorbates/name.adsorabte

        Adsorbate file should read e.g.
        
        ```
        Adsorbate: CH2CH2
        bead,x,y,z,charge
        CH2,0,0,0,0
        CH2,1.33,0,0,0
        ```
        """
        adsorbate = new()
        
        # read adsorbate file
        if ~ isfile("data/adsorbates/" * name * ".adsorbate")
            error(@sprintf("Adsorbate file %s.adsorbate not present in data/adsorbates/", name))
        end
        f = open("data/adsorbates/" * name * ".adsorbate")
        lines = readlines(f)

        adsorbate.nbeads = length(lines) - 2
        adsorbate.name = split(lines[1])[2]
        
        adsorbate.bead_names = AbstractString[]
        adsorbate.bead_xyz = zeros(Float64, 3, adsorbate.nbeads)
        adsorbate.bead_charges = zeros(Float64, adsorbate.nbeads)
        adsorbate.charged_flag = false # change later if true

        assert(lines[2] == "bead,x,y,z,charge\n")

        for i = 1:adsorbate.nbeads
            line = replace(lines[i+2], "\n", "")
            line = split(line, ",")

            push!(adsorbate.bead_names, line[1])

            adsorbate.bead_xyz[1, i] = parse(Float64, line[2])
            adsorbate.bead_xyz[2, i] = parse(Float64, line[3])
            adsorbate.bead_xyz[3, i] = parse(Float64, line[4])

            adsorbate.bead_charges[i] = parse(Float64, line[5])
            if abs(adsorbate.bead_charges[i]) > 0.0000001
                adsorbate.charged_flag = true
            end
        end

        close(f)
    
        # Check that the adsorbate has no net charge.
        if abs(sum(adsorbate.bead_charges)) > .0001
            error(@sprintf("Adsorbate has net charge of %f\n", sum(adsorbate.bead_charges)))
        end

        adsorbate.get_MW = function()
            """
            Get molecular weight of adsorbate.
            Reads molecular mases from data/atomicmasses.csv
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
            Get center of mass (COM) of adsorbate molecule
            Reads molecular mases from data/atomicmasses.csv
            """
            if ! isfile("data/atomicmasses.csv")
                error("Atomic masses file data/atomicmasses.csv not present")
            end
            df_masses = readtable("data/atomicmasses.csv")
            
            mass = 0.0  # total mass
            x_COM = zeros(Float64, 3)
            for i = 1:adsorbate.nbeads
                if ~ (adsorbate.bead_names[i] in df_masses[:atom])
                    error(@sprintf("Adsorbate bead %s not present in data/atomicmasses.cv file", 
                                adsorbate.bead_names[i]))
                end
                
                mass_bead = df_masses[df_masses[:atom] .== adsorbate.bead_names[i], :][:mass][1]
                mass += mass_bead
                x_COM += mass_bead * adsorbate.bead_xyz[:, i]
            end
            return x_COM / mass
        end
        
        # Get center of mass and store it
        adsorbate.x_center_of_mass = adsorbate._get_COM()

        adsorbate.print_bond_lengths = function()
            """
            Print all bond lengths
            """
            @printf("Bond lengths for %s:\n", adsorbate.name)
            @printf("       ")
            for i = 1:adsorbate.nbeads
                @printf("%s     ", adsorbate.bead_names[i])
            end
            for i = 1:adsorbate.nbeads
                @printf("\n%s     ", adsorbate.bead_names[i])
                for j = 1:adsorbate.nbeads
                    if i == j
                        @printf("---- ")
                        continue
                    end
                    bond_length = norm(adsorbate.bead_xyz[:, i] - adsorbate.bead_xyz[:, j])
                    @printf("%.4f ", bond_length)
                end
                @printf("\n")
            end
        end

        adsorbate.center_at_origin = function()
            """
            Set bead_xyz at the center of mass
            """
            for i = 1:adsorbate.nbeads
                adsorbate.bead_xyz[:, i] = adsorbate.bead_xyz[:, i] - adsorbate.x_center_of_mass
            end
            
            # update center of mass
            adsorbate.x_center_of_mass = [0.0, 0.0, 0.0]

            @test_approx_eq_eps adsorbate._get_COM()[1] 0.0 1e-6
            @test_approx_eq_eps adsorbate._get_COM()[2] 0.0 1e-6
            @test_approx_eq_eps adsorbate._get_COM()[3] 0.0 1e-6
        end
        
        # Center adsorbate so its center of mass is at zero
        adsorbate.center_at_origin()
        
        adsorbate.translate_to = function(x::Array{Float64})
            """
            Translate adsorbate in its given configuration from origin to Cartesian vector x
            i.e shift adsorbate so that its center of mass is at x
            """
            for i = 1:adsorbate.nbeads
                # subtract COM to put at origin, then add x
                adsorbate.bead_xyz[:, i] = adsorbate.bead_xyz[:, i] - adsorbate.x_center_of_mass + x
            end
            # update center of mass
            adsorbate.x_center_of_mass = deepcopy(x)
        end 
        
        adsorbate.translate_by = function(dx::Array{Float64})
            """
            Translate adsorbate BY cartesian vector dx
            i.e shift adsorbate so that its center of mass is shifted by dx
            """
            for i = 1:adsorbate.nbeads
                # subtract COM to put at origin, then add x
                adsorbate.bead_xyz[:, i] += dx
            end
            # update center of mass
            adsorbate.x_center_of_mass += dx
        end 
        
        adsorbate.perform_uniform_random_rotation = function ()
            """
            Rotate the adsorbate into a random rotation.
              Precisely, pick an atom. The location of this atom will be moved to a position uniformly 
              distributed on a sphere about the center of mass.
            Keep current center of mass
            # see http://www.mech.utah.edu/~brannon/public/rotation.pdf pg 106
            """
            # build rotation matrix R
            R = zeros(3, 3)
            R[:, 1] = uniform_unit_vector_on_sphere()
            m = uniform_unit_vector_on_sphere()
            R[:, 2] = m - dot(m, R[:, 1]) * R[:, 1]  # subtract of component along col1 so it is orthogonal
            R[:, 2] = R[:, 2] / norm(R[:, 2])
            R[:, 3] = cross(R[:, 1], R[:, 2])  # gives orthogonal vector to first two

            for i = 1:adsorbate.nbeads
                # Note that we need to center adsorbate before multiplying by a rotation matrix
                adsorbate.bead_xyz[:, i] = R * (adsorbate.bead_xyz[:, i] - adsorbate.x_center_of_mass) + adsorbate.x_center_of_mass
            end
        end

        adsorbate.rotate_about_axis = function(theta::Float64, axis::Int)
            """
            Rotate adsorbate molecule (in its given configuration) about a given axis
            Parameters:
                theta: angle of rotation, in degrees
                axis: 1, 2, or 3 for rotating about x-, y-, or z-axis, respectively.
            """
            R = rotation_matrix(theta, axis)
            for i = 1:adsorbate.nbeads
                adsorbate.bead_xyz[:, i] = R * (adsorbate.bead_xyz[:, i] - adsorbate.x_center_of_mass) + adsorbate.x_center_of_mass
            end
        end

        adsorbate.write_to_xyz = function(filename::AbstractString)
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

        return adsorbate
    end
end

function print_info(adsorbate::Adsorbate)
    @printf("Adsorbate name: %s\n", adsorbate.name)
    for i = 1:adsorbate.nbeads
        @printf("Bead %d. x=(%f, %f, %f) charge = %f\n", i, adsorbate.bead_xyz[1, i], adsorbate.bead_xyz[2, i], adsorbate.bead_xyz[3, i], adsorbate.bead_charges[i])
    end
    @printf("center of mass = (%f, %f, %f)", adsorbate.x_center_of_mass[1], adsorbate.x_center_of_mass[2], adsorbate.x_center_of_mass[3])
end
