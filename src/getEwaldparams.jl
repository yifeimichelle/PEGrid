include("framework.jl")

function getEwaldparams(framework::Framework, sr_cutoff::Float64, precision::Float64; verboseflag::Bool=false)
    """
    Goal is to get the following EWald summation parameters:
        alpha: convergence parameter (1/A^2)
        k_reps: number of replications in k-space
    
    from the following:
        framework: crystal structure
        sr_cutoff: short-range cutoff (A)
        precision: EWald precision desired
    """
    ###
    ###  Rules to determine EWald parameters for a given precision (from DLPoly)
    ###  http://www.ccp5.ac.uk/DL_POLY_CLASSIC/FAQ/FAQ2.shtml
    ###
    blah = sqrt(abs(log(precision * sr_cutoff)))
    # alpha parameter to determine width of Gaussian
    alpha = sqrt(abs(log(precision * sr_cutoff * blah))) / sr_cutoff
    # number of replications in k-space for long-range interactions
    blah_ = sqrt(-log(precision * sr_cutoff * (2.0 * blah * alpha)^2))
    k_reps = Array(Int64, 3)
    k_reps[1] = round(Int64, 0.25 + framework.a * alpha * blah_ / pi)
    k_reps[2] = round(Int64, 0.25 + framework.b * alpha * blah_ / pi)
    k_reps[3] = round(Int64, 0.25 + framework.c * alpha * blah_ / pi)
    alpha = alpha^2 # Our alpha is different
    if verboseflag
        @printf("alpha convergence parameter = %f, k_reps = [%d, %d, %d] for Ewald precision %f and sr_cutoff = %f A\n", 
                alpha, k_reps[1], k_reps[2], k_reps[3], precision, sr_cutoff)
    end
    return Dict("alpha" => alpha, 
                "k_reps" => k_reps)
end
