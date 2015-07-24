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
    
    # atom identites
    atoms::Array{String}
    
    # fractional coordinates (3 by natoms)
    fractional_coords::Array{Float64}

    # charges on atoms (units: electron charge)
    charges::Array{Float64}

    # transformation matrix from fractional to cartesian
    f_to_cartesian_mtrx::Array{Float64}
    cartesian_to_f_mtrx::Array{Float64}

    # reciprocal lattice vectors
    reciprocal_lattice::Array{Float64}

    # functions
    crystaldensity::Function
    chemicalformula::Function
    check_for_atom_overlap::Function
    check_for_charge_neutrality::Function
    print_info::Function

    # constructor
    function Framework(structurename::String)
        """
        :param String structurename: name of structure. Will try to import file structurename.cssr
        """
        framework = new()
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
        v_unit_piped = sqrt(1.0 - cos(framework.alpha)^2 - cos(framework.beta)^2 - cos(framework.gamma)^2 +
                2 * cos(framework.alpha) * cos(framework.beta) * cos(framework.gamma)) # volume of unit parallelpiped
        framework.v_unitcell = v_unit_piped * framework.a * framework.b * framework.c  # volume of unit cell
        framework.f_to_cartesian_mtrx = Array(Float64, (3,3))
        framework.cartesian_to_f_mtrx = Array(Float64, (3,3))


        framework.f_to_cartesian_mtrx[1,1] = framework.a
        framework.f_to_cartesian_mtrx[1,2] = framework.b * cos(framework.gamma)
        framework.f_to_cartesian_mtrx[1,3] = framework.c * cos(framework.beta)
        framework.f_to_cartesian_mtrx[2,1] = 0.0
        framework.f_to_cartesian_mtrx[2,2] = framework.b * sin(framework.gamma)
        framework.f_to_cartesian_mtrx[2,3] = framework.c * (cos(framework.alpha) - 
                                             cos(framework.beta) * cos(framework.gamma)) / sin(framework.gamma)
        framework.f_to_cartesian_mtrx[3,1] = 0.0
        framework.f_to_cartesian_mtrx[3,2] = 0.0
        framework.f_to_cartesian_mtrx[3,3] = framework.c * v_unit_piped / sin(framework.gamma)

        framework.cartesian_to_f_mtrx[1,1] = 1.0 / framework.a;
        framework.cartesian_to_f_mtrx[1,2] = - cos(framework.gamma) / framework.a / sin(framework.gamma);
        framework.cartesian_to_f_mtrx[1,3] = (cos(framework.alpha) * cos(framework.gamma) - cos(framework.beta)) / 
                                             (framework.a * v_unit_piped * sin(framework.gamma));
        framework.cartesian_to_f_mtrx[2,1] = 0.0;
        framework.cartesian_to_f_mtrx[2,2] = 1.0 / framework.b / sin(framework.gamma);
        framework.cartesian_to_f_mtrx[2,3] = (cos(framework.beta) * cos(framework.gamma) - cos(framework.alpha)) / 
                                             (framework.b * v_unit_piped * sin(framework.gamma));
        framework.cartesian_to_f_mtrx[3,1] = 0.0;
        framework.cartesian_to_f_mtrx[3,2] = 0.0;
        framework.cartesian_to_f_mtrx[3,3] = sin(framework.gamma) / (framework.c * v_unit_piped);

        # get reciprocal lattice vectors
        framework.reciprocal_lattice = Array(Float64, (3,3))
        framework.reciprocal_lattice[:, 1] = 2.0 * pi * cross(framework.f_to_cartesian_mtrx[:, 2], framework.f_to_cartesian_mtrx[:, 3]) / 
                dot(framework.f_to_cartesian_mtrx[:, 1], cross(framework.f_to_cartesian_mtrx[:, 2], framework.f_to_cartesian_mtrx[:, 3]))
        framework.reciprocal_lattice[:, 2] = 2.0 * pi * cross(framework.f_to_cartesian_mtrx[:, 3], framework.f_to_cartesian_mtrx[:, 1]) / 
                dot(framework.f_to_cartesian_mtrx[:, 2], cross(framework.f_to_cartesian_mtrx[:, 3], framework.f_to_cartesian_mtrx[:, 1]))
        framework.reciprocal_lattice[:, 3] = 2.0 * pi * cross(framework.f_to_cartesian_mtrx[:, 1], framework.f_to_cartesian_mtrx[:, 2]) / 
                dot(framework.f_to_cartesian_mtrx[:, 3], cross(framework.f_to_cartesian_mtrx[:, 1], framework.f_to_cartesian_mtrx[:, 2]))
        
        # get atom count, initialize arrays holding coords
        framework.natoms = int(split(readline(f))[1])
        framework.atoms = Array(String, framework.natoms)
        framework.charges = Array(Float64, framework.natoms)
        framework.fractional_coords = zeros(Float64, 3, framework.natoms)  # fractional coordinates

        # read in atoms and fractional coordinates
        readline(f) # waste a line
        for a = 1:framework.natoms
            line = split(readline(f))

            framework.atoms[a] = line[2]

            framework.fractional_coords[1, a] = float(line[3]) % 1.0 # wrap to [0,1]
            framework.fractional_coords[2, a] = float(line[4]) % 1.0
            framework.fractional_coords[3, a] = float(line[5]) % 1.0

            framework.charges[a] = float(line[14])
        end
        
        close(f) # close file
        
        framework.crystaldensity = function ()
            """
            Computes crystal density of the framework (kg/m3)
            """
            if ! isfile("data/atomicmasses.csv")
                print("Atomic masses file data/atomicmasses.csv not present")
            end
            df = readtable("data/atomicmasses.csv")

            mass = 0.0 # count mass of all atoms here
            for a = 1:framework.natoms
                if ~ (framework.atoms[a] in df[:atom])
                    error(@sprintf("Framework atom %s not present in data/atomicmasses.cv file", 
                                framework.atoms[a]))
                end
                mass += df[df[:atom] .== framework.atoms[a], :][:mass][1]
            end
            
            return mass / framework.v_unitcell * 1660.53892  # --> kg/m3
        end  # end crystaldensity

        framework.chemicalformula = function ()
            """
            Get chemical formula of structure
            """
            # use dictionary to count atom types
            atom_dict = Dict(unique(framework.atoms), zeros(Int, length(unique(framework.atoms))))
            for i = 1:framework.natoms
                atom_dict[framework.atoms[i]] += 1
            end

            # get greatest common divisor
            gcd_ = gcd([k for k in values(atom_dict)]...)
            
            # turn into chemical formula
            for a in keys(atom_dict)
                atom_dict[a] = atom_dict[a] / gcd_
            end

            # print result 
            @printf("Chemical formula of %s:\n\t", framework.structurename)
            for a in keys(atom_dict)
                @printf("%s_%d", a, atom_dict[a])
            end
            @printf("\n")

            # write to file
            if ! isdir(homedir() * "/PEGrid_output/chemicalformulas")
               mkdir(homedir() * "/PEGrid_output/chemicalformulas") 
            end
            formulafile = open(homedir() * "/PEGrid_output/chemicalformulas/" * framework.structurename * ".formula", "w")
            @printf(formulafile, "Atom,Number\n")
            for a in keys(atom_dict)
                @printf(formulafile, "%s,%d\n", a, atom_dict[a])
            end
            close(formulafile)

            return atom_dict
        end  # end chemicalformula

        framework.check_for_atom_overlap = function (distance_tol::Float64)
            """
            Loop through all atoms, check that they do not overlap with any others
            i.e. check for duplicate atoms in the crystal structure file
            :param Float64 distance_tol: tolerance for when an atom is considered to overlap
            :returns false if no overlap
            """
            for i = 1:framework.natoms
                x_i = framework.f_to_cartesian_mtrx * framework.fractional_coords[:, i]
                # loop over adjacent unit cells
                for i_x = -1:1
                    for i_y = -1:1
                        for i_z = -1:1
                            for j = 1:framework.natoms
                                if (j == i)
                                    continue
                                end
                                x_j = framework.f_to_cartesian_mtrx * (framework.fractional_coords[:, j] + 1.0 * [i_x, i_y, i_z])
                                r = norm(x_i - x_j)
                                if (r < distance_tol)
                                    @printf("WARNING: overlap found")
                                    @printf("Atom %d and %d are a distance %f apart.\n", i, j, r)
                                    return true
                                end
                            end  # loop over atom j
                        end  # loop over z uc
                    end  # loop over y uc
                end  # loop over x uc
            end  # loop over atom i
            return false
        end  # end check_for_atom_overlap

        framework.check_for_charge_neutrality = function(tol::Float64)
            """
            Check for charge neutrality within a tolerance
            """
            net_charge = sum(framework.charges)
            if (net_charge < -abs(tol)) | (net_charge > abs(tol))
                @printf("Warning: Framework is not charge neutral!\n")
                @printf("Net charge is %f.\n", net_charge)
            end
        end

        framework.print_info = function()
            @printf("Framework: %s\n", framework.structurename)
            @printf("\t%d atoms.\n", framework.natoms)
            for i = 1:framework.natoms
                @printf("%d. %s, xf = (%f, %f, %f), charge = %f\n", i, framework.atoms[i], 
                        framework.fractional_coords[1, i], framework.fractional_coords[2, i], framework.fractional_coords[3, i],
                        framework.charges[i])
            end
            @printf("Crystal density: %f kg/m3\n", framework.crystaldensity())
            @printf("Total charge: %f\n", sum(framework.charges))
        end

        framework.check_for_atom_overlap(0.1)
        framework.check_for_charge_neutrality(0.001)

        return framework
    end  # end constructor
end  # end Framework type

function shift_unit_cell_write_cssr(frameworkname::String, xf_shift::Array{Float64})
    """
    Shift unit cell by fractional amount and write to a new .cssr
    Shift by fractional amount xf_shift = [.25, .25, .25] for example
    Saves to data/structures/frameworkname_shifted.cssr
    """
    framework = Framework(frameworkname)
    @assert(sum(xf_shift .< 0.0) == 0)
    @assert(sum(xf_shift .> 1.0) == 0)

    f = open("data/structures/" * framework.structurename * "_shifted.cssr", "w")
    write(f, @sprintf("%f %f %f\n", framework.a, framework.b, framework.c))
    write(f, @sprintf("%f %f %f\n", framework.alpha * 180 / pi, framework.beta * 180 / pi, framework.gamma * 180 / pi))
    write(f, @sprintf("%d 0\n", framework.natoms))
    write(f, @sprintf("0 %s : %s\n", framework.structurename, framework.structurename))
    for i = 1:framework.natoms
        # search in all directions
        for i_x = -1:1
            for i_y = -1:1
                for i_z = -1:1
                    # only store if this shifted atom is in [0,1]^3
                    pos = framework.fractional_coords[:, i] + xf_shift + [i_x, i_y, i_z]
                    if ((sum(pos .> 1.0) == 0) & (sum(pos .< 0.0) == 0))
                        write(f, @sprintf("%d %s %f %f %f  0  0  0  0  0  0  0  0  0.000000\n", i, framework.atoms[i], pos[1], pos[2], pos[3]))
                    end
                end
            end
        end
    end
    @printf("Shifted unit cell present in /data/structures/%s_shifted.cssr\n", framework.structurename)
    close(f)
end
    
function write_unitcell_boundary_vtk(frameworkname::String)
    """
    Write unit cell boundary as a .vtk file for visualizing the unit cell.
    """
    framework = Framework(frameworkname)
    if ! isdir(homedir() * "/PEGrid_output")
       mkdir(homedir() * "/PEGrid_output") 
    end
    vtk_file = open(homedir() * "/PEGrid_output/" * framework.structurename * ".vtk", "w")
    
    # write first lines
    @printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\n
                       ASCII\nDATASET POLYDATA\nPOINTS 8 double\n") 

    # write points on boundary of unit cell
    for i = 0:1
        for j = 0:1
            for k = 0:1
                x = [i, j, k] # fractional coordinate
                cellpoint = framework.f_to_cartesian_mtrx * x
                @printf(vtk_file, "%.3f %.3f %.3f\n", cellpoint[1], cellpoint[2], cellpoint[3])
            end
        end
    end

    # define connections
    @printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n
                       2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")
    close(vtk_file)
    @printf(".vtk available at: %s\n", homedir() * "/PEGrid_output/" * framework.structurename * ".vtk")
end

function replicate_cssr_to_xyz(frameworkname::String; rep_factors::Array{Int} = [1, 1, 1])
    """
    Converts a .cssr crystal structure file to .xyz
    Replicates unit cell into rep_factor by rep_factor by rep_factor supercell.
    The primitive unit cell will be in the middle.

    :param String frameworkname: name of crystal structure
    :param Array{Int} rep_factor: number of times to replicate unit cell
    """
    framework = Framework(frameworkname)
    if ! isdir(homedir() * "/PEGrid_output")
       mkdir(homedir() * "/PEGrid_output") 
    end
    xyz_file = open(homedir() * "/PEGrid_output/" * framework.structurename * ".xyz", "w")
    n = framework.natoms * (2 * rep_factors[1] + 1) * (2 * rep_factors[2] + 1) * (2 * rep_factors[3] + 1)
    @printf(xyz_file, "%d\n\n", n)
    
    for a = 1:framework.natoms
        # print atoms in the core unit cell
        xyz = framework.f_to_cartesian_mtrx * framework.fractional_coords[:, a]
        @printf(xyz_file, "%s %f %f %f\n", framework.atoms[a], 
                xyz[1], xyz[2], xyz[3])
    end

    # print atoms in neighboring unit cells
    for a = 1:framework.natoms
        for rep_x = -rep_factors[1]:rep_factors[1]
            for rep_y = -rep_factors[2]:rep_factors[2]
                for rep_z = -rep_factors[3]:rep_factors[3]
                    if (rep_x == 0) & (rep_y == 0) & (rep_z == 0)
                        continue
                    end
                    # fractional coord
                    x_f = framework.fractional_coords[:, a] + 1.0 * [rep_x, rep_y, rep_z]  
                    xyz = framework.f_to_cartesian_mtrx * x_f
                    @printf(xyz_file, "%s %f %f %f\n", framework.atoms[a], 
                            xyz[1], xyz[2], xyz[3])
                end
            end
        end
    end
    close(xyz_file)
    @printf("Replicated .xyz in %s\n", homedir() * "/PEGrid_output/" * framework.structurename * ".xyz")
end
