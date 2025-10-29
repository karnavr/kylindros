include("../functions.jl")
println("Avaliable threads: ", Threads.nthreads())

## FILE PATHS

# ONE SOLUTION
solutions, constants, metadata = readSolution("/home/karnav/Documents/kylindros/results/library/fu/fu_2410222138_N36.jld2");

# VARYING n
results_path = "/home/karnav/Documents/kylindros/results/experiments/fu_vary_coeffs_tol1e13/"

savepath = "plotting scripts/results"

## SINGLE BRANCH PLOTS

# profiles
plot_profiles(solutions, constants, lw = 2.5)
savefig(joinpath(savepath, "fuProfiles.pdf"))

# branch
plot_branch(solutions, constants, color_by_error = false)
savefig(joinpath(savepath, "fuBranch.pdf"))

# profiles with coeffs 
plot_profilesVScoeffs(solutions, constants, indices = [2, 100, 200], shift_profiles = true, n_periods = 2)
savefig(joinpath(savepath, "fuProfilesVSCoeffs.pdf"))

# SINGLE BRANCH 2

solutions2 = "/home/karnav/Documents/kylindros/results/experiments/fu_grid_search/fu_2410004652_N36.jld2"

solutions2, constants2, metadata2 = readSolution(solutions2)

# profiles
plot_profiles(solutions2, constants2, lw = 2.5)
savefig(joinpath(savepath, "fuProfiles2.pdf"))

# branch
plot_branch(solutions2, constants2, color_by_error = false)
savefig(joinpath(savepath, "fuBranch2.pdf"))

## AGGREGATE PLOTS

# cauchy error
plot_cauchy_error(results_path, indices = [50, 100, 150, 175, 200], yscale = :log10)
savefig(joinpath(savepath, "fuCauchyError.pdf"))

# condition number
plot_metric_vs_N(results_path, indices = [1, 50, 75, 100, 125, 150, 175, 200], metric = :condition, legend_position = :topleft)
savefig(joinpath(savepath, "fuConditionNumber.pdf"))

# a collection of different branches for different parameter values

## OTHER