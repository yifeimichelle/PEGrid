PEGrid
======

This code, written in Julia, is for visualizing the potential energy contours of an adsorbate molecule inside a nanoporous material. We superimpose onto the unit cell of the crystal a 3D grid of points and compute the potential energy of an adsorbate molecule at each point. The output of the program is a Gaussian cube file (.cube) of the potential energies (kJ/mol) of the adsorbate at the grid points.

With this .cube energy grid file, we can visualize the potential energy contours of the adsorbate inside the pores of the crystal as in the figure below.

<a href="url"><img src="https://www.dropbox.com/s/2kdo64m8yq092e9/example.png?dl=1" align="center" height="500" width="500" ></a>

## Necessary data

The data folder contains the crystal structure file, the force field used to compute the potential energy of the adsorbate, and the atomic masses of the [pseudo-]atoms.

#### Crystal structure file (.cssr)

Put the crystal structure files in 

#### Force field file (.csv)

#### Atomic masses file (.csv)

## How to write the grid

As an example, if we want to use the `UFF` forcefield to compute the energy of adsorbate molecule `CH4` (modeled as a Lennard-Jones sphere) in crystal structure `IRMOF-1` using a Lennard-Jones cutoff of 12.5 Angstrom on a 3D grid of points with a spacing of 1.0 Angstrom, the following two lines of code will write a .cube grid file (units: kJ/mol) to your home directory.

    include("src/energygrid.jl")
    writegrid("CH4", "IRMOF-1", "UFF", gridspacing=1.0,cutoff=12.5)

* `CH4`: corresponds to an atom listed in data/forcefields/UFF.csv
* `IRMOF-1`: corresponds to crystal structure file data/structures/IRMOF-1.cssr

## Other features

#### Replicating a .cssr to an .xyz file

#### Computing the density of a crystal structure

#### Generating a .vtk mesh file of edges of unit cell 
