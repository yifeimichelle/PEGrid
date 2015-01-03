PEGrid
======

This code, written in Julia, is for visualizing the potential energy contours of an adsorbate molecule inside a nanoporous material. We superimpose onto the unit cell of the crystal a 3D grid of points and compute the potential energy of an adsorbate molecule at each point. The output of the program is a Gaussian cube file (file extension: .cube) of the potential energies (units: kJ/mol) of the adsorbate at the grid points.

With this .cube energy grid file, we can visualize the potential energy contours of the adsorbate inside the pores of the crystal as in the figure below.

<a href="url"><img src="https://www.dropbox.com/s/2kdo64m8yq092e9/example.png?dl=1" align="center" height="500" width="500" ></a>

## Necessary data

The data folder contains the crystal structure file, the force field used to compute the potential energy of the adsorbate, and the atomic masses of the [pseudo-]atoms.

PEGrid requires the `DataFrames` package in Julia. To install, type in Julia: `Pkg.add("DataFrames")`.

#### Crystal structure file (.cssr)

PEGrid reads crystal structure files in the .cssr format. If needed to convert .cif to .cssr, use Zeo++. 

    ./network -cssr ${yourstructurename}.cif

Place the crystal structure files in `data/structures`.

#### Force field file (.csv)

The force field is the model and parameters used to describe the potential energy of the adsorbate molecule with the atoms of the crystal structure. PEGrid models the interaction between the adsorbate *a* and crystal structure atom type *i* a distance *r* apart with the Lennard-Jones potential:

![](https://www.dropbox.com/s/zfn2titfiyp8w6j/LJpotential.png?dl=1)

Currently, only adsorbates modeled as a Lennard-Jones sphere are supported; no electrostatic charges or more complex molecules.

The Lennard-Jones parameters for the cross-interaction between adsorbate *a* and atom type *i* are computed from the pure *a-a*, *i-i* interactions using the following Lorenz-Berthelot mixing rules:

![](https://www.dropbox.com/s/zfn2titfiyp8w6j/LJpotential.png?dl=1)

The pure *i-i* and *a-a* Lennard-Jones interaction parameters must be stored in `data/forcefields/${yourforcefield}.csv`.

We mimic an infinite crystal by applying periodic boundary conditions. This is enabled by approximating interactions beyond a cutoff radius to be zero. The cutoff radius, technically part of the force field, is provided as an argument in the functions of PEGrid.

#### Atomic masses file (.csv)

In order to calculate the crystal density of the framework, PEGrid stores the atomic masses (units: amu) of the atoms in `data/atomicmasses.csv`.

## How to write the grid

As an example, if we want to use the `UFF` forcefield to compute the energy of adsorbate molecule `CH4` (modeled as a Lennard-Jones sphere) in crystal structure `IRMOF-1` using a Lennard-Jones cutoff of 12.5 Angstrom on a 3D grid of points with a spacing of 1.0 Angstrom, the following two lines of code in Julia will write a .cube grid file (units: kJ/mol) to your home directory.

    include("src/energygrid.jl")
    writegrid("CH4", "IRMOF-1", "UFF", gridspacing=1.0, cutoff=12.5)

* `CH4`: corresponds to an atom listed in data/forcefields/UFF.csv
* `IRMOF-1`: corresponds to crystal structure file data/structures/IRMOF-1.cssr

I recommend using the IJulia notebook. The code will print off the progress of the grid writing every 10%. Be patient, as computing fine grids and/or large unit cells lead to long computation times.

## How to visualize potential energy contours with the grid

Use VisIt visualization tool.

- [ ] Interface with PEGrid

## Other features

#### Replicating a .cssr to an .xyz file

An .xyz file contains a list of atoms in the crystal structure file, their identities, and their positions in Cartesian coordinates. An .xyz file is useful for visualizing the atoms of the crystal structure file. PEGrid can replicate the unit cell of the crystal in each direction, so that the 'home' unit cell is in the center, and output a .xyz of the replicated crystal structure file.

The function `replicate_cssr_to_xyz` will replicate the unit cell of the crystal and save an .xyz file in the current directory (just above the `data` folder).

    include("src/framework.jl")
    replicate_cssr_to_xyz(${yourstructurename})

#### Computing the density of a crystal structure

To compute the crystal density of a framework, first create a `Framework` type, which reads the crystal coordinates, atomic identities, and unit cell sizes from the crystal structure file:

    include("src/framework.jl")
    framework = constructframework("${yourstructurename}")

Then, call the function `crystaldensity(framework::Framework)`

    rho = crystaldensity(framework)

#### Generating a .vtk mesh file of edges of unit cell 

The .vtk mesh file allows us to visualize the boundaries of the unit cell. To generate a .vtk of the unit cell boundary, use the following Julia code:

    include("src/framework.jl")
    write_unitcell_boundary_vtk("${yourstructurename}")

A `${yourstructurename}.vtk` will be written in the working directory (just above the `data` directory).

#### Computing the potential energy of an adsorbate at a particular point

e.g., to compute the potential energy of adsorbate `CH4` in crystal structure `IRMOF-1` using Lennard-Jones parameters in `data/forcefield/UFF.csv` with a cutoff radius of 12.5 at the *fractional* coordinate `(x_f, y_f, z_f)`, use the Julia code:

    include("src/energyutils.jl")
    E = E_vdw_at_point("IRMOF-1", "UFF", "CH4", x_f, y_f, z_f, cutoff=12.5)

- [ ] include function to convert from Cartesian to fractional for this function
