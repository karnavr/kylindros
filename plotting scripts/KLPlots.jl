include("../functions.jl")

## FILE PATHS

# ONE SOLUTION
solutions, constants, metadata = readSolution("/home/karnav/Documents/kylindros/results/library/kirchhoffLove/kirchhoffLove_2710223338_N36.jld2");

# VARYING n
results_path = "/home/karnav/Documents/kylindros/results/experiments/kl_vary_coeffs/"

savepath = "plotting scripts/results"

## SINGLE BRANCH PLOTS

# profiles
plot_profiles(solutions, constants, lw = 2.5)
savefig(joinpath(savepath, "klProfiles.pdf"))

# branch
plot_branch(solutions, constants, color_by_error = false)
savefig(joinpath(savepath, "klBranch.pdf"))

# profiles with coeffs 
plot_profilesVScoeffs(solutions, constants, indices = [2, 50, 100], shift_profiles = true, n_periods = 2)
savefig(joinpath(savepath, "klProfilesVSCoeffs.pdf"))


## AGGREGATE PLOTS

# cauchy error
plot_cauchy_error(results_path, indices = [50, 100, 150, 175, 200], yscale = :log10)
savefig(joinpath(savepath, "klCauchyError.pdf"))

# condition number
plot_metric_vs_N(results_path, indices = [50, 75, 100, 125, 150, 175, 200], metric = :condition, legend_position = :topleft)
savefig(joinpath(savepath, "klConditionNumber.pdf"))

# a collection of different branches for different parameter values

## OTHER