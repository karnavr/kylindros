# Parameter grid search for EB model varying L

# Load the project-wide functions and dependencies
include(joinpath(@__DIR__, "..", "functions.jl"))

using CSV, DataFrames

# Fixed parameters
branchN = 200
a1Vals = collect(range(0.01, 0.23, branchN + 1))
N = 36
b = 0.1

# Parameter sweep ranges
L_values = collect(range(π/2, 8π, 20))

save_dir = "experiments/eb_grid_search"

# Create error log file
error_log_path = joinpath(dirname(@__DIR__), "results", "experiments", "eb_grid_search", "error_log.csv")
mkpath(dirname(error_log_path))

# Initialize error log DataFrame
error_log = DataFrame(L = Float64[], error_type = String[], error_message = String[], timestamp = String[])

@info "Starting parameter grid search" total_runs=length(L_values) date=Dates.now()

completed = 0
failed = 0

for (i, L) in enumerate(L_values)
    @info "Run $i/$(length(L_values))" L=L
    
    try
        # Create constants
        constants = eulerBernoulliConstants(N, L, b)

        # Initialize wave speed from k1
        k1 = 1 * π / constants.L
        cInitial = c0(k1, constants)

        # Initial guess
        initial_guess = (1e-16) .* ones(branchN + 1, constants.N + 2)
        initial_guess[1, 1:4] = [cInitial, 1.0, a1Vals[1], 1e-16]

        # Compute branch
        @time begin
            solutions, constants_out, metadata = bifurcation(
                initial_guess,
                a1Vals,
                branchN,
                constants,
                tol = 1e-14,
                solver = :NLSolver,
                max_iter = 1000,
                overwrite = true,
                save_dir = save_dir,
                verbose = false,
            )
            # Explicitly drop references
            solutions = nothing
            constants_out = nothing
            metadata = nothing
        end

        global completed += 1
        @info "Completed run" L=L
        
    catch err
        global failed += 1
        error_type = string(typeof(err))
        error_message = replace(sprint(showerror, err), '\n' => ' ')
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        
        # Append to error log
        push!(error_log, (L, error_type, error_message, timestamp))
        CSV.write(error_log_path, error_log)
        flush(stdout)  # Ensure logs are written immediately
        
        @warn "Error during run" L=L error=error_type
        
    finally
        # Proactively free memory between runs
        GC.gc()
    end
end

@info "Grid search complete" completed=completed failed=failed date=Dates.now()

if failed > 0
    @info "Error log saved to" path=error_log_path
end