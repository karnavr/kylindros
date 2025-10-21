# Standalone script to compute and plot Cauchy error from existing fuSimple results
# Usage: plot_cauchy_error(results_dir, indices=[50, 100, 150])

include(joinpath(@__DIR__, "..", "functions.jl"))

"""
    get_profile_on_grid(coeffs, z_grid, L)

Reconstruct surface profile S(z) from Fourier coefficients on arbitrary grid.
"""
function get_profile_on_grid(coeffs, z_grid, L)
    return fourierSeries(coeffs, z_grid, L)[1]
end

"""
    map_files_by_N(dirpath)

Walk directory and return Dict mapping N => filepath for each available N.
"""
function map_files_by_N(dirpath)
    if !isdir(dirpath)
        error("Directory not found: $dirpath")
    end
    
    files_by_N = Dict{Int, String}()
    
    for (root, _, filenames) in walkdir(dirpath)
        for fname in filenames
            if endswith(fname, ".jld2")
                fpath = joinpath(root, fname)
                try
                    _, constants, _ = readSolution(fpath)
                    N = Int(constants.N)
                    # Keep latest file for each N (or just overwrite)
                    files_by_N[N] = fpath
                catch
                    # Skip unreadable files
                    continue
                end
            end
        end
    end
    
    return files_by_N
end

"""
    compute_cauchy_error(fileN_path, fileNm2_path, idx)

Compute Cauchy error E_N = ||S_N - S_{N-2}||_2 for branch index idx.
Both profiles are evaluated on the finer z_N grid.
"""
function compute_cauchy_error(fileN_path, fileNm2_path, idx)
    # Load both solutions
    solsN, constsN, metaN = readSolution(fileN_path)
    solsNm2, constsNm2, metaNm2 = readSolution(fileNm2_path)
    
    # Sanity checks
    if constsN.L != constsNm2.L
        error("L mismatch between N=$(constsN.N) and N=$(constsNm2.N)")
    end
    
    if idx > size(solsN, 1) || idx > size(solsNm2, 1)
        error("Index $idx out of bounds")
    end
    
    # Extract coefficients
    coeffsN = solsN[idx, 2:end]
    coeffsNm2 = solsNm2[idx, 2:end]
    
    # Use finer grid from N run
    z_fine = constsN.z
    L = constsN.L
    
    # Reconstruct both profiles on same grid
    S_N = get_profile_on_grid(coeffsN, z_fine, L)
    S_Nm2 = get_profile_on_grid(coeffsNm2, z_fine, L)
    
    # Compute L2 norm
    diff_sq = (S_N .- S_Nm2).^2
    E_N = sqrt(trapz(z_fine, diff_sq))
    
    return E_N
end

"""
    plot_cauchy_error(results_dir; indices=[50, 100, 150], N_range=nothing)

Main function: plot Cauchy error vs N for specified branch indices.

# Arguments
- `results_dir`: Path to results directory (e.g., "results/fuSimpleVaryCoeffs/fuSimple")
- `indices`: Vector of branch indices to plot (each gets a separate curve)
- `N_range`: Optional range of N values to include (default: all available)

# Returns
Plots object showing log10(E_N) vs N for each branch index.
"""
function plot_cauchy_error(results_dir; indices=[50, 100, 150], N_range=nothing)
    
    # Map all available files by N
    files_by_N = map_files_by_N(results_dir)
    
    if isempty(files_by_N)
        error("No valid result files found in $results_dir")
    end
    
    # Get sorted list of available N values
    N_vals = sort(collect(keys(files_by_N)))
    
    # Filter by N_range if provided
    if N_range !== nothing
        N_vals = filter(n -> n in N_range, N_vals)
    end
    
    println("Found N values: $N_vals")
    
    # Storage for plot data: idx => [(N, E_N, a1), ...]
    data_by_idx = Dict{Int, Vector{Tuple{Int, Float64, Float64}}}()
    for idx in indices
        data_by_idx[idx] = []
    end
    
    # Compute Cauchy errors for each N where N-2 exists
    for N in N_vals
        if N - 2 in keys(files_by_N)
            fileN = files_by_N[N]
            fileNm2 = files_by_N[N - 2]
            
            # Try to compute for each requested index
            for idx in indices
                try
                    E_N = compute_cauchy_error(fileN, fileNm2, idx)
                    
                    # Get a1 value for labeling (from higher-N file)
                    sols, _, _ = readSolution(fileN)
                    a1 = sols[idx, 3]
                    
                    # Store if valid
                    if E_N > 0 && isfinite(E_N)
                        push!(data_by_idx[idx], (N, E_N, a1))
                    end
                catch e
                    @warn "Failed to compute Cauchy error" N=N idx=idx error=e
                    continue
                end
            end
        end
    end
    
    # Create plot
    p = plot(legend=:topright, size=(700, 500))
    
    # Plot each branch index as a separate series
    for idx in indices
        data = data_by_idx[idx]
        
        if isempty(data)
            @warn "No data for index $idx"
            continue
        end
        
        # Sort by N
        perm = sortperm(data; by = x -> x[1])
        sorted = data[perm]
        
        N_plot = [t[1] for t in sorted]
        E_plot = [t[2] for t in sorted]
        a1_label = sorted[end][3]  # Use a1 from largest N
        
        # Plot with log scale
        plot!(N_plot, E_plot,
            yaxis = :log10,
            lw = 2,
            marker = :circle,
            markersize = 4,
            label = "a₁ = $(round(a1_label, digits=3))", 
            size = (500, 500))
    end
    
    xlabel!("N")
    ylabel!(L"E_N")
    
    return p
end

# Example usage (uncomment to run):
# results_path = "/home/karnav/Documents/kylindros/results/fuSimpleVaryCoeffs/fuSimple"
# p = plot_cauchy_error(results_path; indices=[50, 100, 150])
# display(p)

