# read cif from Avogadro

type Atom
    element::AbstractString
    x::Float64
    y::Float64
    z::Float64
end

function read_Avogadro_cif(filename::AbstractString)
    """
    Read .cif file saved from Avogadro.

    The problem is that it stores Cartesian coords instead of fractional
    """
    data = Dict()
    fields = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
              "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]
    
    f = open("data/structures/cifs/" * filename)
    atoms = Atom[]
    readcoords = false
    for line in eachline(f)
        for i = 1:length(fields)
            if split(line)[1] == fields[i]
                data[fields[i]] = parse(Float64, split(line)[2])
            end
            if split(line)[1] == "_space_group_name_Hall"
                println(split(line))
                @assert(split(line)[2:3] == ["'P", "1'"])
            end
        end

        if readcoords == true
            element = split(line)[1]
            x = parse(Float64, split(line)[3])
            y = parse(Float64, split(line)[4])
            z = parse(Float64, split(line)[5])
            push!(atoms, Atom(element, x, y, z))
        end
        
        if split(line)[1] == "_atom_site_Cartn_z"
            readcoords = true
        end
    end
    close(f)
   
    # build transformation matrix 
    a = data["_cell_length_a"]
    b = data["_cell_length_b"]
    c = data["_cell_length_c"]

    alpha = data["_cell_angle_alpha"]
    beta = data["_cell_angle_beta"]
    gamma = data["_cell_angle_gamma"]
        
    v_unit_piped = sqrt(1.0 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 +
                2 * cos(alpha) * cos(beta) * cos(gamma)) # volume of unit parallelpiped

    cartesian_to_f_mtrx = Array(Float64, (3,3))

    cartesian_to_f_mtrx[1,1] = 1.0 / a
    cartesian_to_f_mtrx[1,2] = - cos(gamma) / a / sin(gamma)
    cartesian_to_f_mtrx[1,3] = (cos(alpha) * cos(gamma) - cos(beta)) / 
                                         (a * v_unit_piped * sin(gamma))
    cartesian_to_f_mtrx[2,1] = 0.0
    cartesian_to_f_mtrx[2,2] = 1.0 / b / sin(gamma)
    cartesian_to_f_mtrx[2,3] = (cos(beta) * cos(gamma) - cos(alpha)) / 
                                         (b * v_unit_piped * sin(gamma))
    cartesian_to_f_mtrx[3,1] = 0.0
    cartesian_to_f_mtrx[3,2] = 0.0
    cartesian_to_f_mtrx[3,3] = sin(gamma) / (c * v_unit_piped)
    
    # write new cif with fractional
    f = open("data/structures/cifs/" * replace(filename, ".cif", "") * "_.cif", "w")
    @printf(f, "_cell_length_a %f\n", a)
    @printf(f, "_cell_length_b %f\n", b)
    @printf(f, "_cell_length_c %f\n", c)

    @printf(f, "_cell_angle_alpha %f\n", alpha)
    @printf(f, "_cell_angle_beta %f\n", beta)
    @printf(f, "_cell_angle_gamma %f\n", gamma)
    
    @printf(f, "\n_symmetry_space_group_name_H-M      'P1'\n")
    @printf(f, "_symmetry_Int_Tables_number     1\n")
    @printf(f, "\nloop_\n_symmetry_equiv_pos_as_xyz\n'+x,+y,+z'\n\nloop_\n")
    @printf(f, "_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")

    for i = 1:length(atoms)
        xf = cartesian_to_f_mtrx * [atoms[i].x, atoms[i].y, atoms[i].z]
        @printf(f, "%d %s %f %f %f\n", i, atoms[i].element,
                        xf[1], xf[2], xf[3])
    end

    close(f)
    @printf("New cif with fractional coords present in %s\n", "data/structures/cifs/" * replace(filename, ".cif", "") * "_.cif")

    return atoms
end
