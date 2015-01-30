function evaluate_expressions_on_all_cores(;kwargs...)
    """
    Evaluates a given expression on each process. For broadcasting variables
    """
    for (nm, val) in kwargs # unpacks expressions in args
        for p = 1:nprocs()
            @spawnat p eval(Main, Expr(:(=), nm, val))
        end
    end
end

function tester(a::Int)
    """
    Compute the potential energy of an adsorbate molecule on a 3D grid of points superimposed on the unit cell of the structure.
    Parallelized across cores.

    The grid is written to a file `structurename.cube`, in Gaussian cube format. The units of the energy are kJ/mol.

    :param: String adsorbate: the name of the adsorbate molecule, corresponding to the forcefield file
    """
    @everywhere function geta()
        a
    end

    # broadcast function arguments to all cores
    evaluate_expressions_on_all_cores(a=a)
    for p in workers()
        r = @spawnat p geta()
        @printf("value stored on process %d: %d\n", p, fetch(r))
    end

end

tester(3)
