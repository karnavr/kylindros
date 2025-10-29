# Batch compute KL branches for N = 8:48 with robust error handling

# Load the project-wide functions and dependencies
include(joinpath(@__DIR__, "..", "functions.jl"))

# Parameters that are constant across runs
branchN = 200
a1Vals = collect(range(0.01, 0.29, branchN + 1))

L = π
b = 0.1

save_dir = "experiments/kl_vary_coeffs"

@info "Starting batch run for N=8:48" date=Dates.now()

for N in 8:48
    @info "Starting N run" N=N
    try
        # Create constants
        constants = kirchhoffLoveConstants(N, L, b)

        # Initialize wave speed from k1
        k1 = 1 * π / constants.L
        cInitial = c0(k1, constants)

        # Initial guess
        initial_guess = (1e-16) .* ones(branchN + 1, constants.N + 2)
        initial_guess[1, 1:4] = [cInitial, 1.0, a1Vals[1], 1e-16]

        # Compute branch (overwrite to ensure deterministic re-runs)
        @time begin
            solutions, constants_out, metadata = bifurcation(
                initial_guess,
                a1Vals,
                branchN,
                constants,
                tol = 1e-13,
                solver = :NLSolver,
                max_iter = 1000,
                overwrite = true,
                save_dir = save_dir,
            )
            # Explicitly drop references to large results immediately after save
            solutions = nothing
            constants_out = nothing
            metadata = nothing
        end

        @info "Completed N run" N=N
    catch err
        # Log and continue with next N
        @warn "Error during N run, skipping" N=N error=err
    finally
        # Proactively free memory between runs
        GC.gc()
    end
end

@info "Batch run complete" date=Dates.now()

