include("framework.jl")
using Roots

function getEwaldparams(framework::Framework, precision::Float64)
    """
    alpha: convergence parameter (1/A^2)
    rc: short range cutoff (A)
    k_reps: number of replications in k-space

    From Understanding Molecular Simulation by Frenkel and Smit
    Chapter 12.1 Long range interactions
    """
    sqrt_alpha = (framework.natoms * pi^3 / framework.v_unitcell^2) ^ (1 / 6)
    # k I think there is a mistake in Berend's book... 
    # this alpha is different from one in Ewald derivation by square

    g(x) = exp(- x^2) / x^2 - precision
    s = fzero(g, 0.1)
    rc = s / sqrt_alpha

    k_reps = [framework.a, framework.b, framework.c] * s * sqrt_alpha / pi

    @printf("rc = %f A,  alpha = %f 1/A^2, kreps = [%f, %f, %f]\n", rc, sqrt_alpha ^ 2, k_reps[1], k_reps[2], k_reps[3])
    return sqrt_alpha ^ 2, rc, k_reps
end
