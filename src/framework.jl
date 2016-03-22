###
#   Author: Cory M. Simon (CoryMSimon@gmail.com)
###
 # module Framework
 # 
 # export Framework, write_to_cssr, replicate_cssr_to_xyz, shift_unit_cell, consolidate_atoms_with_same_charge, write_unitcell_boundary_vtk, replicate_cssr_to_xyz

using DataFrames


function read_crystal_file(structurename::AbstractString, file_extension::AbstractString)
    """
    Read crystal structure file.
    file_extension in ["cif", "cssr"]
    """
    @assert(file_extension in ["cif", "cssr"], "Reads only cif or cssr")

    data = Dict()  # store all data here and return later

    atoms = AbstractString[]
    xf = Float64[]
    yf = Float64[]
    zf = Float64[]
    charges = Float64[]
    
    # open crystal structure file and read its lines
    f = open("data/structures/" * structurename * "." * file_extension)
    if ~ isfile("data/structures/" * structurename * "." * file_extension)
        @printf("Crystal structure file %s not present in data/structures/", structurename * "." * file_extension)
    end
    lines = readlines(f)

    if file_extension == "cif"
        loop_starts = -1
        for i = 1:length(lines)
            line = split(lines[i])
            if length(line) == 0
                continue
            end
            if line[1] == "_symmetry_space_group_name_H-M"
                @assert(contains(line[2] * line[3], "P1"), ".cif must have P1 symmetry.\n")
            end
            
            for axis in ["a", "b", "c"]
                if line[1] == @sprintf("_cell_length_%s", axis)
                    data[axis] = parse(Float64, line[2])
                end
            end
            for angle in ["alpha", "beta", "gamma"]
                if line[1] == @sprintf("_cell_angle_%s", angle)
                    data[angle] = parse(Float64, line[2]) * pi / 180.0
                end
            end

            if (line[1] == "loop_")
                next_line = split(lines[i+1])
                if (next_line[1] == "_atom_site_label")
                    loop_starts = i + 1
                    break
                end
            end
        end  # end loop over lines

        # broke the loop. so loop_starts is line where "_loop" first starts
        name_to_column = Dict{AbstractString, Int}()
        i = loop_starts
        while length(split(lines[i])) == 1
            name_to_column[split(lines[i])[1]] = i + 1 - loop_starts
            i += 1
        end
        println(name_to_column)

        # now extract fractional coords of atoms and their charges
        for i = loop_starts+length(name_to_column):length(lines)
            line = split(lines[i])
            push!(atoms, line[name_to_column["_atom_site_label"]])
            push!(xf, mod(parse(Float64, line[name_to_column["_atom_site_fract_x"]]), 1.0))
            push!(yf, mod(parse(Float64, line[name_to_column["_atom_site_fract_y"]]), 1.0))
            push!(zf, mod(parse(Float64, line[name_to_column["_atom_site_fract_z"]]), 1.0))
            push!(charges, parse(Float64, line[name_to_column["_atom_site_charge"]]))
            if length(line) != length(name_to_column)
                break
            end
        end
        data["natoms"] = length(xf)
    end  # if file is .cif

    if file_extension == "cssr"
        # get unit cell sizes
        line = split(lines[1]) # first line is (a, b, c)
        data["a"] = parse(Float64, line[1])
        data["b"] = parse(Float64, line[2])
        data["c"] = parse(Float64, line[3])

        # get unit cell angles. Convert to radians
        line = split(lines[2])
        data["alpha"] = parse(Float64, line[1]) * pi / 180.0
        data["beta"] = parse(Float64, line[2]) * pi / 180.0
        data["gamma"] = parse(Float64, line[3]) * pi / 180.0

        data["natoms"] = parse(Int, split(lines[3])[1])

        # read in atoms and fractional coordinates
        for a = 1:data["natoms"]
            line = split(lines[4+a])
            
            push!(atoms, line[2])

            push!(xf, mod(parse(Float64, line[3]), 1.0)) # wrap to [0,1]
            push!(yf, mod(parse(Float64, line[4]), 1.0)) # wrap to [0,1]
            push!(zf, mod(parse(Float64, line[5]), 1.0)) # wrap to [0,1]

            push!(charges, parse(Float64, line[14]))
        end
    end

    close(f) # close file

    data["charges"] = charges
    data["atoms"] = atoms
    data["fractional_coords"] = transpose(hcat(xf, yf, zf))
    return data
end

type Framework
    """
    Stores info of a crystal structure

    Constructed by reading a .cssr
    """
    structurename::AbstractString

    "unit cell Bravais lattice size (Angstrom)"
    a::Float64
    b::Float64
    c::Float64

    # unit cell Bravais lattice angles
    alpha::Float64
    beta::Float64
    gamma::Float64

    # volume of the unit cell (Angstrom^3)
    v_unitcell::Float64

    # atom count (in unit cell)
    natoms::Int
    
    # atom identites
    atoms::Array{AbstractString}
    
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
    m_unitcell::Function  # mass of the unit cell (amu)
    chemicalformula::Function
    check_for_atom_overlap::Function
    check_for_charge_neutrality::Function
    print_info::Function
    get_COM::Function
    write_to_cssr::Function
    write_to_cif::Function
    replicate_framework_to_xyz::Function
    write_unitcell_boundary_vtk::Function
    append_number_labels_to_atoms::Function
    remove_number_labels_on_atoms::Function

    # constructor
    function Framework(structurename::AbstractString, check_for_atom_overlap::Bool=true; file_extension::AbstractString="cssr")
        """
        :param AbstractString structurename: name of structure. Will try to import file structurename.cssr
        """
        framework = new()
        framework.structurename = structurename
        
        data = read_crystal_file(structurename, file_extension)

        # get unit cell sizes
        framework.a = data["a"]
        framework.b = data["b"]
        framework.c = data["c"]

        # get unit cell angles (Radians)
        framework.alpha = data["alpha"]
        framework.beta = data["beta"]
        framework.gamma = data["gamma"]

        # transformation matrices for fractional coords <--> cartesian coords
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
        framework.natoms = data["natoms"]
        framework.atoms = data["atoms"]
        framework.charges = data["charges"]
        framework.fractional_coords = data["fractional_coords"]
        @assert(size(framework.fractional_coords) == (3, framework.natoms))

        framework.m_unitcell = function()
            """
            Get mass of unit cell (amu)
            """
            # get atomic mass dictionary
            if ! isfile("data/atomicmasses.csv")
                print("Atomic masses file data/atomicmasses.csv not present")
            end
            df = readtable("data/atomicmasses.csv")
            mass_dict = Dict()
            for i = 1:size(df, 1)
                mass_dict[df[:atom][i]] = df[:mass][i] 
            end
            
            # get mass of unit cell
            mass = 0.0 # count mass of all atoms here
            for a = 1:framework.natoms
                if ~ (haskey(mass_dict, framework.atoms[a]))
                    error(@sprintf("Framework atom %s not present in data/atomicmasses.cv file", 
                                framework.atoms[a]))
                end
                mass += mass_dict[framework.atoms[a]]
            end
            return mass  # amu
        end
        
        framework.crystaldensity = function ()
            """
            Computes crystal density of the framework (kg/m3)
            """
            return framework.m_unitcell() / framework.v_unitcell * 1660.53892  # --> kg/m3
        end  # end crystaldensity

        framework.chemicalformula = function ()
            """
            Get chemical formula of structure
            """
            # use dictionary to count atom types
            atom_dict = Dict(zip(unique(framework.atoms), zeros(Int, length(unique(framework.atoms)))))
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
            overlap_flag = false
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
                                    overlap_flag = true
                                end
                            end  # loop over atom j
                        end  # loop over z uc
                    end  # loop over y uc
                end  # loop over x uc
            end  # loop over atom i
            return overlap_flag
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

        framework.get_COM = function()
            """
            Compute center of mass
            """
            com = [0.0, 0.0, 0.0]
            # create mass dictionary
            if ! isfile("data/atomicmasses.csv")
                print("Atomic masses file data/atomicmasses.csv not present")
            end
            df = readtable("data/atomicmasses.csv")
            mass_dict = Dict()
            for i = 1:size(df, 1)
                mass_dict[df[:atom][i]] = df[:mass][i] 
            end

            # get masses
            masses = map(x -> mass_dict[x], framework.atoms)

            # wrap fractional coords onto circle
            eta = cos(framework.fractional_coords * 2 * pi)
            zeta = sin(framework.fractional_coords * 2 * pi)
            # take mean
            eta_hat = mean(eta, 2)
            zeta_hat = mean(zeta, 2)
            tol = 1e-5  # tolerance for declaring zero (atan unstable)
            for i = 1:3
                if ((eta_hat[i] < tol) & (eta_hat[i] > -tol))
                    eta_hat[i] = 0.0
                end
                if ((zeta_hat[i] < tol) & (zeta_hat[i] > -tol))
                    zeta_hat[i] = 0.0
                end
                com[i] = (atan2(-zeta_hat[i], -eta_hat[i]) + pi) / (2.0 * pi)
            end

            return com
        end

        framework.write_to_cssr = function(filename::AbstractString)
            """
            Write to .cssr with desired filename
            """
            if (filename == framework.structurename * ".cssr")
                error("With this filename, we will overwrite the original structure...\n")
            end
            f = open("data/structures/" * filename, "w")
            write(f, @sprintf("\t\t\t%f %f %f\n", framework.a, framework.b, framework.c))
            write(f, @sprintf("\t\t%f %f %f SPGR = 1 P 1      OPT = 0\n", 
                                framework.alpha * 180 / pi, 
                                framework.beta * 180 / pi, 
                                framework.gamma * 180 / pi)
            )
            write(f, @sprintf("%d 0\n", framework.natoms))
            write(f, @sprintf("0 %s : %s\n", framework.structurename, "revised by PEGrid"))
            for i = 1:framework.natoms
                # only store if this shifted atom is in [0,1]^3
                write(f, @sprintf(" %d %s %f %f %f  0  0  0  0  0  0  0  0  %f\n", i, framework.atoms[i], 
                        framework.fractional_coords[1, i], framework.fractional_coords[2, i], framework.fractional_coords[3, i],
                        framework.charges[i]))
            end
            @printf("Cssr can be found in /data/structures/%s\n", filename)
            close(f)
        end

        framework.write_to_cif = function(filename::AbstractString)
            """
            Write framework to .cif format
            """
            if (filename == framework.structurename * ".cif")
                error("With this filename, we will overwrite the original structure...\n")
            end
            f = open("data/structures/" * filename * ".cif", "w")
            @printf(f, "_symmetry_space_group_name_H-M   'P 1'\n")

            @printf(f, "_cell_length_a %f\n", framework.a)
            @printf(f, "_cell_length_b %f\n", framework.b)
            @printf(f, "_cell_length_c %f\n", framework.c)
            
            @printf(f, "_cell_angle_alpha %f\n", framework.alpha * 180.0 / pi)
            @printf(f, "_cell_angle_beta %f\n", framework.beta * 180.0 / pi)
            @printf(f, "_cell_angle_gamma %f\n", framework.gamma * 180.0 / pi)
            @printf(f, "_cell_volume %f\n", framework.v_unitcell)

            @printf(f, "_symmetry_Int_Tables_number 1\n\n")
            @printf(f, "loop_\n_symmetry_equiv_pos_as_xyz\n 'x, y, z'\n\n")

            @printf(f, "loop_\n_atom_site_label\n_atom_site_type_symbol\n")
            @printf(f, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
            @printf(f, "_atom_site_charge\n")

             for i = 1:framework.natoms
                # elements may be e.g Ca21
                element = framework.atoms[i]
                if ! isalpha(element)
                    if isalpha(element[1:2])
                        element = element[1:2]
                    else
                        element = element[1]
                    end
                    @assert(isalpha(element))
                end
                @printf(f, "%s %s %f %f %f %f\n", framework.atoms[i], element,
                            framework.fractional_coords[1, i],
                            framework.fractional_coords[2, i],
                            framework.fractional_coords[3, i],
                            framework.charges[i])
             end
             close(f)
             @printf("%s.cif file present in data/structures/\n", filename)
        end
        
        framework.write_unitcell_boundary_vtk = function()
            """
            Write unit cell boundary as a .vtk file for visualizing the unit cell.
            """
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

        framework.replicate_framework_to_xyz = function(rep_factors::Array{Int}, alldirections::Bool)
            """
            Converts a Framework to an xyz file; replications of unit cell tunable.
            Replicates unit cell into rep_factor by rep_factor by rep_factor supercell.
            The home unit cell will be in the middle.

            :param Array{Int} rep_factor: number of times to replicate unit cell
            :param Bool alldirections: replicate in all directions if true. else, just 0 to positive reps
            """
            if ! isdir(homedir() * "/PEGrid_output")
               mkdir(homedir() * "/PEGrid_output") 
            end
            xyz_file = open(homedir() * "/PEGrid_output/" * framework.structurename * ".xyz", "w")
            if alldirections
                n = framework.natoms * (2 * rep_factors[1] + 1) * (2 * rep_factors[2] + 1) * (2 * rep_factors[3] + 1)
                lower_bound_reps = -rep_factors
            else
                n = framework.natoms * (rep_factors[1] + 1) * (rep_factors[2] + 1) * (rep_factors[3] + 1)
                lower_bound_reps = [0, 0, 0]
            end
            @printf(xyz_file, "%d\n\n", n)
            
            for a = 1:framework.natoms
                # print atoms in the core unit cell
                xyz = framework.f_to_cartesian_mtrx * framework.fractional_coords[:, a]
                @printf(xyz_file, "%s %f %f %f\n", framework.atoms[a], 
                        xyz[1], xyz[2], xyz[3])
            end

            # print atoms in neighboring unit cells
            for a = 1:framework.natoms
                for rep_x = lower_bound_reps[1]:rep_factors[1]
                    for rep_y = lower_bound_reps[2]:rep_factors[2]
                        for rep_z = lower_bound_reps[3]:rep_factors[3]
                            if (rep_x == 0) & (rep_y == 0) & (rep_z == 0)
                                # already printed this above
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

        framework.append_number_labels_to_atoms = function()
            """
            Appends numbers to each atom.
            e.g. C->C1, C->C2
            """
            atoms = unique(f.atoms)
            for atom in atoms
                idx = find(f.atoms .== atom)
                for i = 1:length(idx)
                    f.atoms[idx[i]] = @sprintf("%s%d", atom, i)
                end
            end
        end

        framework.remove_number_labels_on_atoms = function()
            """
            Removes number labels from atoms
            e.g. C1->C, C2->C
            """
            for i = 1:framework.natoms
                while !isalpha(f.atoms[i][end])
                    f.atoms[i] = chop(f.atoms[i])
                end
            end
        end

        if check_for_atom_overlap
            framework.check_for_atom_overlap(0.1)
        end
        framework.check_for_charge_neutrality(0.001)

        return framework
    end  # end constructor
end  # end Framework type

function shift_unit_cell(framework::Framework, xf_shift::Array{Float64})
    """
    Shift unit cell by fractional amount
    Shift by fractional amount xf_shift = [.25, .25, .25] for example
    """
    new = deepcopy(framework)
    @assert(sum(xf_shift .< 0.0) == 0)
    @assert(sum(xf_shift .> 1.0) == 0)
    
    count = 0
    for i = 1:framework.natoms
        # search in all directions
        for i_x = -1:1
            for i_y = -1:1
                for i_z = -1:1
                    # only store if this shifted atom is in [0,1]^3
                    pos = framework.fractional_coords[:, i] + xf_shift + 1.0 * [i_x, i_y, i_z]
                    if ((sum(pos .> 1.0) == 0) && (sum(pos .< 0.0) == 0))
                        count += 1
                        println(count)
                        new.atoms[count] = framework.atoms[i]
                        new.fractional_coords[1, count] = pos[1]
                        new.fractional_coords[2, count] = pos[2]
                        new.fractional_coords[3, count] = pos[3]
                    end
                end
            end
        end
    end
    @assert(count == framework.natoms)
    return new
end

function consolidate_atoms_with_same_charge(framework::Framework; decimal_tol::Int=3, verbose::Bool=false)
    """
    Groups atoms of same element and charge together as a single type.
    e.g. C with charge = 0.5 ==> C_a
    decimal_tol = 2 considers 0.231 and 0.234 the same charge
    returns new_framework, charge_dictionary, element_dictionary
    """
    new = deepcopy(framework)
    
    strings_to_append = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"]
    # get unique atom types
    unique_atoms = unique(framework.atoms)
    # get rounded charges for comparison
    rounded_charges = round(framework.charges, decimal_tol)
    if verbose
        @printf("New_atom_label Element Charge\n")
    end

    # charge dictionary. atom_label:charge
    charge_dict = Dict()
    # element dictionary. atom_label:element
    element_dict = Dict()
    # for every atom, make it C_a etc
    for a in unique_atoms
        idx = framework.atoms .== a
        # get unique charges of this atom
        unique_charges = unique(rounded_charges[idx])
        # consolidate atoms a with the same charge
        for c = 1:length(unique_charges)
            if c > length(strings_to_append)
                error("need more strings to append to element names...")
            end
            new_label = a * "_" * strings_to_append[c]
            
            idx_ = idx & (rounded_charges .== unique_charges[c])
            new.atoms[idx_] = new_label
            if verbose
                @printf("%s %s %s\n", new_label, a, unique_charges[c])
            end
            element_dict[new_label] = a
            charge_dict[new_label] = unique_charges[c]
        end
    end
    new.charges = rounded_charges
    
    # ensure charge neutrality
    shift_charges_by = sum(new.charges) / new.natoms
    new.charges = new.charges - shift_charges_by
    
    # update charges
    for elem in keys(charge_dict)
        charge_dict[elem] -= shift_charges_by
    end

    return new, charge_dict, element_dict
end
    
