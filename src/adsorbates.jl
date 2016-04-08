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

    # Store coordinates here, column-wise
    x::Array{Float64}  # Cartesian
    base_x::Array{Float64}  # base Cartesian coords

    bead_names::Array{AbstractString}  # corresponds to name in force field
    charges::Array{Float64}  # charges on each bead
    x_center_of_mass::Array{Float64}  # store center of mass so we don't need to compute it each time

    translate_to::Function  # translate adsorbate such that its center of mass is at x
    translate_to_from_base::Function  # translate adsorbate such that its center of mass is at x
    translate_by::Function # translate adsorbate by dx
    perform_uniform_random_rotation::Function  # perform a rotation of the adsorbate
    write_to_xyz::Function  # write adsorbate positions to .xyz
    rotate_about_axis::Function  # rotate adsorbate about a given axis
    translate_and_rotate_to::Function  # rotate given Euler angles and position (about base coords)
    print_bond_lengths::Function

    molecular_weight::Function
    _get_COM::Function
    center_at_origin ::Function  # subtract center of mass so that it is at origin

    function Adsorbate(name::AbstractString)
        """
        Constructor. Information for adsorbate is stored in file:
            PEGRID_DATA_DIR/adsorbates/name.adsorabte

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
        if ~ isfile(PEGRID_DATA_DIR * "/adsorbates/" * name * ".adsorbate")
            error(@sprintf("Adsorbate file %s.adsorbate not present in %s/adsorbates/", name, PEGRID_DATA_DIR))
        end
        f = open(PEGRID_DATA_DIR * "/adsorbates/" * name * ".adsorbate")
        lines = readlines(f)

        adsorbate.nbeads = length(lines) - 2
        adsorbate.name = split(lines[1])[2]
        
        adsorbate.bead_names = AbstractString[]
        adsorbate.x = zeros(Float64, 3, adsorbate.nbeads)
        adsorbate.base_x = zeros(Float64, 3, adsorbate.nbeads)
        adsorbate.charges = zeros(Float64, adsorbate.nbeads)
        adsorbate.charged_flag = false # change later if true

        assert(lines[2] == "bead,x,y,z,charge\n")

        for i = 1:adsorbate.nbeads
            line = replace(lines[i+2], "\n", "")
            line = split(line, ",")

            push!(adsorbate.bead_names, line[1])
            
            adsorbate.x[1, i] = parse(Float64, line[2])
            adsorbate.x[2, i] = parse(Float64, line[3])
            adsorbate.x[3, i] = parse(Float64, line[4])

            adsorbate.charges[i] = parse(Float64, line[5])
            if abs(adsorbate.charges[i]) > 0.0000001
                adsorbate.charged_flag = true
            end
        end

        close(f)
    
        # Check that the adsorbate has no net charge.
        if abs(sum(adsorbate.charges)) > .0001
            error(@sprintf("Adsorbate has net charge of %f\n", sum(adsorbate.charges)))
        end

        adsorbate.molecular_weight = function()
            """
            Get molecular weight of adsorbate.
            Reads molecular mases from data/atomicmasses.csv
            """
            if ! isfile(PEGRID_DATA_DIR * "/atomicmasses.csv")
                error(@sprintf("Atomic masses file %s/atomicmasses.csv not present.\n", PEGRID_DATA_DIR))
            end
            df_masses = readtable(PEGRID_DATA_DIR * "/atomicmasses.csv")
            
            MW = 0.0
            for i = 1:adsorbate.nbeads
                if ~ (adsorbate.bead_names[i] in df_masses[:atom])
                    error(@sprintf("Adsorbate bead %s not present in %s/atomicmasses.csv\n", 
                                adsorbate.bead_names[i], PEGRID_DATA_DIR))
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
            if ! isfile(PEGRID_DATA_DIR * "/atomicmasses.csv")
                error(@sprintf("Atomic masses file %s/atomicmasses.csv not present.\n", PEGRID_DATA_DIR))
            end
            df_masses = readtable(PEGRID_DATA_DIR * "/atomicmasses.csv")
            
            mass = 0.0  # total mass
            x_COM = zeros(Float64, 3)
            for i = 1:adsorbate.nbeads
                if ~ (adsorbate.bead_names[i] in df_masses[:atom])
                    error(@sprintf("Adsorbate bead %s not present in %s/atomicmasses.csv\n", 
                                adsorbate.bead_names[i], PEGRID_DATA_DIR))
                end
                
                mass_bead = df_masses[df_masses[:atom] .== adsorbate.bead_names[i], :][:mass][1]
                mass += mass_bead
                x_COM += mass_bead * adsorbate.x[:, i]
            end
            return x_COM / mass
        end
        
        # Get center of mass and store it
        adsorbate.x_center_of_mass = adsorbate._get_COM()

        adsorbate.center_at_origin = function()
            """
            Translate adsorbate such that its center of mass is at the origin.
            """
            for i = 1:adsorbate.nbeads
                adsorbate.x[:, i] = adsorbate.x[:, i] - adsorbate.x_center_of_mass
            end
            
            # update center of mass
            adsorbate.x_center_of_mass = [0.0, 0.0, 0.0]

            x_com = adsorbate._get_COM()

            @test_approx_eq_eps x_com[1] 0.0 1e-6
            @test_approx_eq_eps x_com[2] 0.0 1e-6
            @test_approx_eq_eps x_com[3] 0.0 1e-6
            return
        end
        
        # Shift adsorbate so its center of mass is at origin
        adsorbate.center_at_origin()

        # Declare these coords as the "base" coords
        adsorbate.base_x = deepcopy(adsorbate.x)
        
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
                    bond_length = norm(adsorbate.x[:, i] - adsorbate.x[:, j])
                    @printf("%.4f ", bond_length)
                end
                @printf("\n")
            end
            return
        end

        adsorbate.translate_to = function(x::Array{Float64})
            """
            Translate adsorbate in its given configuration from origin to Cartesian vector x
            i.e shift adsorbate so that its center of mass is at x
            """
            for i = 1:adsorbate.nbeads
                # subtract COM to put at origin, then add x
                adsorbate.x[:, i] = adsorbate.x[:, i] - adsorbate.x_center_of_mass + x
            end
            # update center of mass
            adsorbate.x_center_of_mass = deepcopy(x)
            return
        end 
        
        adsorbate.translate_to_from_base = function(x::Array{Float64})
            """
            Translate adsorbate in its given configuration from origin to Cartesian vector x
            i.e shift adsorbate so that its center of mass is at x
            Use base coordinates.
            """
            for i = 1:adsorbate.nbeads
                # subtract COM to put at origin, then add x
                adsorbate.x[:, i] = adsorbate.base_x[:, i] + x
            end
            # update center of mass
            adsorbate.x_center_of_mass = deepcopy(x)
            return
        end 
        
        adsorbate.translate_by = function(dx::Array{Float64})
            """
            Translate adsorbate BY cartesian vector dx
            i.e shift adsorbate so that its center of mass is shifted by dx
            """
            for i = 1:adsorbate.nbeads
                # subtract COM to put at origin, then add x
                adsorbate.x[:, i] += dx
            end
            # update center of mass
            adsorbate.x_center_of_mass += dx
            return
        end 
        
        adsorbate.perform_uniform_random_rotation = function ()
            """
            Rotate the adsorbate into a random rotation.
              Precisely, pick an atom. The location of this atom will be moved to a position uniformly 
              distributed on a sphere about the center of mass.
            Keep current center of mass
            # see http://www.mech.utah.edu/~brannon/public/rotation.pdf pg 106
            Since it doesn't matter, we use the base_coords here to avoid numerical drift
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
                # base coords hv COM at origin
                adsorbate.x[:, i] = R * adsorbate.base_x[:, i] + adsorbate.x_center_of_mass
            end
            return
        end

        adsorbate.rotate_about_axis = function(theta::Float64, axis::Int)
            """
            Rotate adsorbate molecule (in its current configuration) about a given axis
            Parameters:
                theta: angle of rotation, in degrees
                axis: 1, 2, or 3 for rotating about x-, y-, or z-axis, respectively.
            """
            R = rotation_matrix(theta, axis)
            for i = 1:adsorbate.nbeads
                adsorbate.x[:, i] = R * (adsorbate.x[:, i] - adsorbate.x_center_of_mass) + adsorbate.x_center_of_mass
            end
            return
        end

        adsorbate.translate_and_rotate_to = function(x::Array{Float64}, alpha::Float64, beta::Float64, gamma::Float64)
            """
            Put COM of adsorbate at x and rotate adsorbate using Euler angles alpha, beta, and gamma. Do this WRT base_x.

            Parameters:
                x: 3D vector to which we translate the adsorbate
                alpha, beta, gamma: Euler angles, in radians
            """
            R_x = rotation_matrix(alpha, 1)
            R_y = rotation_matrix(alpha, 2)
            R_z = rotation_matrix(alpha, 3)

            for i = 1:adsorbate.nbeads
                adsorbate.x[:, i] = R_z * R_y * R_x * adsorbate.base_x[:, i] + x
            end
            return
        end 

        adsorbate.write_to_xyz = function(filename::AbstractString)
            """
            Write adsorbate to .xyz file.

            Atom "C_co2" will be written as "C".
            """
            filename = homedir() * "/PEGrid_output/adsorbates/" * filename * ".xyz"
            
            if ! isdir(homedir() * "/PEGrid_output/adsorbates")
                mkdir(homedir() * "/PEGrid_output/adsorbates")
            end
            f = open(filename, "w")
            write(f, @sprintf("%d\n\n", adsorbate.nbeads))
            for i = 1:adsorbate.nbeads
                atom_write_name = split(adsorbate.bead_names[i], "_")[1]
                write(f, @sprintf("%s %f %f %f\n", atom_write_name, adsorbate.x[1, i], adsorbate.x[2, i], adsorbate.x[3, i]))
            end
            close(f)
            @printf("See %s\n", filename)
            return
        end

        return adsorbate
    end
end

function print_info(adsorbate::Adsorbate)
    @printf("Adsorbate name: %s\n", adsorbate.name)
    for i = 1:adsorbate.nbeads
        @printf("Bead %d. x=(%f, %f, %f) charge = %f\n", i, adsorbate.x[1, i], adsorbate.x[2, i], adsorbate.x[3, i], adsorbate.charges[i])
    end
    @printf("center of mass = (%f, %f, %f)", adsorbate.x_center_of_mass[1], adsorbate.x_center_of_mass[2], adsorbate.x_center_of_mass[3])
end
