"Plotting functions for the results of the numerical continuation."

## Functions for each type of plot

function plot_profiles(solutions, constants; shift_profiles = true, line_color = theme_palette(:default)[1], lw = 2)

	## plot three profiles from the solutions branch (equally spaced in the last 0.75 of the branch)

	# get needed constants
	L = constants.L
	z = constants.z
	branchN = length(solutions[:,1])

	# seperate coeffs and speeds
	coeffs = solutions[:,2:end]

	# pick three profiles with equally spaced a₁ amplitudes from the last 0.75 of the branch
	start_idx = Int(round(0.25*branchN))
	a1_vals = solutions[start_idx:end, 3]  # a₁ values from last 75%
	indices_range = start_idx:branchN
	
	# target three evenly spaced a₁ values
	a1_min = minimum(a1_vals)
	a1_max = maximum(a1_vals)
	target_a1 = range(a1_min, a1_max, length=3)
	
	# find closest actual solution for each target a₁
	indices = [indices_range[argmin(abs.(a1_vals .- target))] for target in target_a1]

	# create array for profiles
	profiles = zeros(branchN,length(z))
	
	# convert profiles
	Threads.@threads for i = 1:branchN
		profiles[i,:] .= fourierSeries(coeffs[i,:], z, L)[1]
	end

	# shift profiles to the right over by L
	if shift_profiles == true
		profiles = [profiles[:,Int(end/2)+1:end] profiles[:,1:Int(end/2)]]; nothing
	end

	linestyles = [:dashdot, :dash, :solid]

	# plot profiles
	p = plot(legend = true, size = (500,500))
	for (i, index) in enumerate(indices)
		plot!(z, profiles[index,:], label = "a₁ = $(round(solutions[index,3], digits=3))", lw=lw, linestyle = linestyles[i], color = line_color)
	end
	xlabel!(L"z"); ylabel!(L"S")

	return p
end

function repeat_profiles(solutions, constants; indices = [size(solutions,1)], n_periods::Int = 3, figure_size = nothing, lw = 2)

	## plot up to three profiles (by branch indices) repeated n_periods times with shared y-axis scale

	# get needed constants
	L = constants.L
	z = constants.z

	# input validation
	if n_periods < 1
		error("n_periods must be ≥ 1")
	end
	if length(indices) < 1 || length(indices) > 3
		error("provide 1 to 3 indices")
	end
	branchN = size(solutions, 1)
	for idx in indices
		if idx < 1 || idx > branchN
			error("index $idx is out of bounds (1:$branchN)")
		end
	end

	# helper to build a repeated profile and extended domain without duplicate junctions
	function build_tiled_profile(coeffs, z, L, n_periods)
		profile = fourierSeries(coeffs, z, L)[1]
		period = z[end] - z[1]
		z_ext = Vector{eltype(z)}()
		prof_ext = Vector{eltype(profile)}()
		for i = 0:(n_periods-1)
			if i == 0
				append!(z_ext, z)
				append!(prof_ext, profile)
			else
				append!(z_ext, (z[2:end] .+ i * period))
				append!(prof_ext, profile[2:end])
			end
		end
		return z_ext, prof_ext
	end

	# build tiled series for requested indices
	z_series = Vector{Vector{eltype(z)}}()
	prof_series = Vector{Vector{Float64}}()
	for idx in indices
		coeffs = solutions[idx, 2:end]
		z_ext, prof_ext = build_tiled_profile(coeffs, z, L, n_periods)
		push!(z_series, z_ext)
		push!(prof_series, prof_ext)
	end

	# compute shared y-limits
	global_min = minimum(minimum.(prof_series))
	global_max = maximum(maximum.(prof_series))
	y_margin = 0.05 * (global_max - global_min)
	if y_margin == 0
		y_margin = 0.01
	end
	ylims_shared = (global_min - y_margin, global_max + y_margin)

	# assemble subplots in a single row
	subplots = []
	for i in 1:length(indices)
		p_i = plot(z_series[i], prof_series[i],
			lw = lw,
			color = :steelblue,
			legend = false,
			ylims = ylims_shared)
		xlabel!(L"z")
		ylabel!(L"S")
		push!(subplots, p_i)
	end

	# stack subplots top-to-bottom and either use provided size or auto-scale height
	if figure_size === nothing
		fig_width = 1000
		fig_height = 300 * length(subplots)
		p = plot(subplots..., layout = (length(subplots), 1), size = (fig_width, fig_height))
	else
		p = plot(subplots..., layout = (length(subplots), 1), size = figure_size)
	end
	return p
end

# Backward-compatible method: single index positional argument
function repeat_profiles(solutions, constants::Constants, index::Int; n_periods::Int = 3, figure_size = nothing)
	return repeat_profiles(solutions, constants; indices = [index], n_periods = n_periods, figure_size = figure_size)
end

function plot_branch(solutions, metadata; color_by_error = false)

	## plot the branch of solutions

	# seperate coeffs and speeds
	speeds = solutions[:,1]
	a₁ = solutions[:,3]

	# plot the branch
	if color_by_error
		errors = metadata["errors"]
		p = scatter(speeds, a₁, 
			marker_z = log10.(abs.(errors)),
			color = :plasma,
			colorbar_title = "log₁₀(Error)",
			colorbar_titlefontsize = 6,
			colorbar_tickfontsize = 6,
			colorbar = false,
			legend = false, 
			size = (500,500))
	else
		p = scatter(speeds, a₁, 
			legend = false, 
			size = (500,500))
	end
	
	xlabel!(L"c")
	ylabel!(L"a_1")

	return p
end

function plot_coeffs(solutions, indices = round.(Int, range(1, size(solutions,1), length=5)))

	## plot the coefficients of the solutions

	# get coeffs
	coeffs = solutions[:,2:end]

	# plot the coefficients
	p = plot(legend=false, size = (500,500))
	for index in indices
		plot!(abs.(coeffs[Int(index),:]), 
			yaxis=:log, 
			label = "a₁ = $(round(solutions[Int(index),3], digits=3))", 
			marker = :circle,
			markersize = 3,
			markerstrokewidth = 0,
			markerstrokealpha = 0)
	end

	# only show legend if there is more than one index
	if length(indices) > 1
		plot!(legend = true)
	end

	xlabel!("n")
	ylabel!(L"a_n")

	return p
end

function plot_error(solutions, metadata)

	# get errors
	errors = metadata["errors"]

	# get a1 values 
	a1 = solutions[:,3]

	# plot error over the branch with both line and markers
	p = plot(a1, errors, 
		lw=2,                   # line width
		marker=:circle,         # marker style
		markersize=2,           # marker size
		yaxis=:log10,           # log scale y-axis
		legend=false)           # no legend

	xlabel!(L"a_1")
	ylabel!("Error")

	return p
end

## Plot everything 

function plotEverything(solutions, constants, metadata)

	branchN = Int(metadata["branchN"])

	profileplot = plot_profiles(solutions, constants, shift_profiles = true)
	branchplot = plot_branch(solutions, metadata)
	coeffplot = plot_coeffs(solutions, 1:(branchN/10):branchN)
	errorplot = plot_error(solutions, metadata)

	p = plot(profileplot, branchplot, coeffplot, errorplot, layout = (2,2), size = (1000,1000))

	return p
end

function plot_comparison_wrapper(solutions, constants::Constants; 
                                indices = [Int(round(0.75*size(solutions,1))), size(solutions,1)],
                                shift_profiles = true,
                                save_figure = false,
                                filename = nothing,
                                figure_size = (1200, 800))

	## Create a figure comparing profiles and coefficients for multiple solutions
	
	# Validate inputs
	if save_figure && filename === nothing
		error("filename must be provided when save_figure=true")
	end
	
	if length(indices) < 1
		error("at least 1 index must be provided")
	end
	
	branchN = size(solutions, 1)
	for idx in indices
		if idx < 1 || idx > branchN
			error("index $idx is out of bounds (1:$branchN)")
		end
	end
	
	# Calculate y-axis limits for consistent scaling across all profiles
	L = constants.L
	z = constants.z
	all_profiles = []
	
	for idx in indices
		coeffs = solutions[idx, 2:end]
		profile = fourierSeries(coeffs, z, L)[1]
		if shift_profiles
			profile = [profile[Int(end/2)+1:end]; profile[1:Int(end/2)]]
		end
		push!(all_profiles, profile)
	end
	
	# Find global y-limits
	global_min = minimum(minimum.(all_profiles))
	global_max = maximum(maximum.(all_profiles))
	y_margin = 0.05 * (global_max - global_min)
	ylims = (global_min - y_margin, global_max + y_margin)
	
	# Helper function to plot a single profile with consistent y-limits
	function plot_single_profile(solutions, constants, index, ylims; shift_profiles = true)
		# get needed constants
		L = constants.L
		z = constants.z
		
		# get coefficients for this solution
		coeffs = solutions[index, 2:end]
		
		# convert to profile
		profile = fourierSeries(coeffs, z, L)[1]
		
		# shift profile if requested
		if shift_profiles
			profile = [profile[Int(end/2)+1:end]; profile[1:Int(end/2)]]
		end
		
		# create plot with consistent y-limits
		p = plot(z, profile, 
			label = "a₁ = $(round(solutions[index,3], digits=3))", 
			lw = 2, 
			color = :steelblue,
			legend = false,
			ylims = ylims,
			size = (400, 400))
		xlabel!(L"z")
		ylabel!(L"S")
		
		return p
	end
	
	# Helper function to plot coefficients for a single solution
	function plot_single_coeffs(solutions, index)
		# get coefficients for this solution
		coeffs = solutions[index, 2:end]
		
		# create plot
		p = plot(abs.(coeffs), 
			yaxis = :log, 
			label = "a₁ = $(round(solutions[index,3], digits=3))", 
			marker = :circle,
			markersize = 3,
			markerstrokewidth = 0,
			markerstrokealpha = 0,
			legend = false,
			size = (400, 400))
		xlabel!("n")
		ylabel!(L"a_n")
		
		return p
	end
	
	# Create plots for all indices
	profile_plots = []
	coeff_plots = []
	
	for idx in indices
		push!(profile_plots, plot_single_profile(solutions, constants, idx, ylims; shift_profiles = shift_profiles))
		push!(coeff_plots, plot_single_coeffs(solutions, idx))
	end
	
	# Strict two-column layout: profiles left, coefficients right
	n_solutions = length(indices)
	all_plots = []
	for i in 1:n_solutions
		push!(all_plots, profile_plots[i])
		push!(all_plots, coeff_plots[i])
	end

	# Create final plot with no empty panes
	p = plot(all_plots..., layout = (n_solutions, 2), size = figure_size)
	
	# Save figure if requested
	if save_figure
		savefig(p, filename)
		println("Figure saved as: $filename")
	end
	
	return p
end


function plot_metric_vs_N(dirpath::String; indices::Vector{Int}, metric::Symbol = :error, yscale::Symbol = :log10, legend_position = :best)

	# Use the directory exactly as provided and gather all .jld2 files recursively
	if !isdir(dirpath)
		error("results directory not found: $(dirpath)")
	end

	files = String[]
	for (root, _, fs) in walkdir(dirpath)
		for f in fs
			endswith(f, ".jld2") && push!(files, joinpath(root, f))
		end
	end

	# Buckets per requested branch index: idx => vector of (N, y, a1)
	data_by_index = Dict{Int, Vector{Tuple{Int, Float64, Float64}}}()
	for idx in indices
		data_by_index[idx] = Vector{Tuple{Int, Float64, Float64}}()
	end

	# Iterate over result files
	for file_path in files
		try
			solutions, constants, metadata = readSolution(file_path)

			branch_len = size(solutions, 1)
			for idx in indices
				if idx <= branch_len
					N = Int(constants.N)
					a1 = solutions[idx, 3]
					if metric == :error
						y = Float64(metadata["errors"][idx])
					elseif metric == :condition
						y = Float64(metadata["condition_numbers"][idx])
					else
						error("unknown metric: $(metric). use :error or :condition")
					end
					# Skip non-positive values on log scale to avoid warnings and empty plots
					if !(yscale == :log10 && y <= 0)
						push!(data_by_index[idx], (N, y, a1))
					end
				end
			end
		catch
			# skip unreadable/broken files silently to keep plotting robust
		end
	end

	# Create plot
	p = plot(legend = legend_position, size = (500, 500))
	for idx in indices
		entries = data_by_index[idx]
		isempty(entries) && continue
		# sort by N
		perm = sortperm(entries; by = x -> x[1])
		sorted = entries[perm]
		xN = [t[1] for t in sorted]
		yvals = [t[2] for t in sorted]
		# use amplitude from largest N as label
		label_a1 = sorted[end][3]
		label = "a₁ = $(round(label_a1, digits=3))"
		plot!(xN, yvals, lw = 2, marker = :circle, markersize = 3, label = label)
	end

	# axes
	xlabel!(L"N")
	if metric == :error
		ylabel!("Error")
	elseif metric == :condition
		ylabel!(L"\sigma (J)")
	end

	if yscale == :log10
		plot!(yaxis = :log10)
	end

	return p
end

function plot_cauchy_error(dirpath::String; indices::Vector{Int}, yscale::Symbol = :log10, legend_position = :best)

	"Plot Cauchy convergence error E_N = ||S_N - S_{N-2}||_{L2} vs N for specified branch indices."

	# ========== Directory validation ==========
	if !isdir(dirpath)
		error("Results directory not found: $(dirpath)")
	end

	# ========== Build file mapping by N ==========
	# Walk directory tree and map N => filepath for each available truncation level
	files_by_N = Dict{Int, String}()

	for (root, _, filenames) in walkdir(dirpath)
		for fname in filenames
			if endswith(fname, ".jld2")
				fpath = joinpath(root, fname)
				try
					_, constants, _ = readSolution(fpath)
					N = Int(constants.N)
					files_by_N[N] = fpath  # Keep latest file for each N
				catch
					# Skip unreadable files silently
					continue
				end
			end
		end
	end

	if isempty(files_by_N)
		error("No valid result files found in $(dirpath)")
	end

	N_vals = sort(collect(keys(files_by_N)))

	# ========== Compute Cauchy errors ==========
	# Storage: idx => [(N, E_N, a1), ...]
	data_by_idx = Dict{Int, Vector{Tuple{Int, Float64, Float64}}}()
	for idx in indices
		data_by_idx[idx] = []
	end

	# Compute errors for each N where N-2 exists
	for N in N_vals
		if (N - 2) in keys(files_by_N)
			fileN = files_by_N[N]
			fileNm2 = files_by_N[N - 2]

			# Try to compute for each requested index
			for idx in indices
				try
					E_N, a1 = compute_cauchy_error(fileN, fileNm2, idx)

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

	# ========== Create plot ==========
	p = plot(legend = legend_position, size = (500, 500))

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

		plot!(N_plot, E_plot,
			lw = 2,
			marker = :circle,
			markersize = 4,
			label = "a₁ = $(round(a1_label, digits=3))")
	end

	# ========== Axes and formatting ==========
	xlabel!(L"N")
	ylabel!(L"E_N")

	if yscale == :log10
		plot!(yaxis = :log10)
	end

	return p
end

## DISPERSION RELATIONS

# ferrofluid
function plot_dispersion(k_range, constants::ferrofluidConstants; vary_param::Symbol = :B, param_range = range(1, 30, length=5))

	"Plots the dispersion relation for a the ferrofluid problem with a varying parameter."

	# unpack constants for the constants struct
	N = constants.N
	L = constants.L
	b = constants.b
	B = constants.B

	# initialize array for speeds
	speeds = zeros(length(param_range), length(k_range))

	for (i, param_value) in enumerate(param_range)

		# create constants struct based on which parameter we are varying
		if vary_param == :B
			constants = ferrofluidConstants(N, L, param_value, b)
		elseif vary_param == :b
			constants = ferrofluidConstants(N, L, B, param_value)
		end
		
		# compute speeds
		# speeds[i,:] .= c0(k_range, constants)

		for (j, k) in enumerate(k_range)
			speeds[i, j] = c0(k, constants)
		end
	end

	# plot speeds
	p = plot(legend=:outertopright)
	param_name = String(vary_param)

	for (i, param_value) in enumerate(param_range)
		plot!(k_range, speeds[i,:], 
			xlabel="Wavenumber (k)", 
			ylabel="Wave speed (c₀)", 
			label="$param_name = $(param_value)", lw=2)
	end

	return p
end

#fuSimple
function plot_dispersion(k_range, constants::fuSimpleConstants; vary_param::Symbol = :λ1, param_range = range(0.1, 2.0, length=5))

	"Plots the dispersion relation for the fu and il'ichev model with a varying parameter."

	# initialize array for speeds
	speeds = zeros(length(param_range), length(k_range))

	for (i, param_value) in enumerate(param_range)

		# create constants struct based on which parameter we are varying
		if vary_param == :λ1
			constants = fuSimpleConstants(constants.N, constants.L, constants.b, param_value, constants.vf)
		elseif vary_param == :λ2
			constants = fuSimpleConstants(constants.N, constants.L, constants.b, constants.λ2, param_value)
		elseif vary_param == :vf
			constants = fuSimpleConstants(constants.N, constants.L, constants.b, constants.λ2, param_value)
		end
		
		# compute speeds
		for (j, k) in enumerate(k_range)
			speeds[i, j] = c0(k, constants)
		end
	end

	# plot speeds
	p = plot(legend=:outertopright)
	param_name = String(vary_param)

	for (i, param_value) in enumerate(param_range)
		plot!(k_range, speeds[i,:], 
			xlabel="Wavenumber (k)", 
			ylabel="Wave speed (c₀)", 
			label="$param_name = $(param_value)", lw=2)
	end

	return p
end 

#fu 
function plot_dispersion(k_range, constants::fuConstants; vary_param::Symbol = :b, param_range = range(0.01, 1.0, length=5))

	"Plots the dispersion relation for the fu and il'ichev model with a varying parameter."

	# initialize array for speeds
	speeds = zeros(length(param_range), length(k_range))

	for (i, param_value) in enumerate(param_range)

		# create constants struct based on which parameter we are varying
		if vary_param == :λ1
			constants = fuConstants(constants.N, constants.L, constants.b, param_value, constants.λ2, constants.vf)
		elseif vary_param == :λ2
			constants = fuConstants(constants.N, constants.L, constants.b, constants.λ1, param_value, constants.vf)
		elseif vary_param == :vf
			constants = fuConstants(constants.N, constants.L, constants.b, constants.λ1, constants.λ2, param_value)
		elseif vary_param == :b
			constants = fuConstants(constants.N, constants.L, param_value, constants.λ1, constants.λ2, constants.vf)
		end
		
		# compute speeds
		# speeds[i,:] .= c0(k_range, constants)

		for (j, k) in enumerate(k_range)
			speeds[i, j] = c0(k, constants)
		end
	end

	# plot speeds
	p = plot(legend=:outertopright)
	param_name = String(vary_param)

	for (i, param_value) in enumerate(param_range)
		plot!(k_range, speeds[i,:], 
			xlabel="Wavenumber (k)", 
			ylabel="Wave speed (c₀)", 
			label="$param_name = $(param_value)", lw=2)
	end

	return p
end

## PARAMETER VARIATION PLOTS

function plotBranchVariations(dir_path::String; 
                             constant_param::Symbol, 
                             constant_value::Real,
                             varied_param_range::Union{Nothing, Tuple{Real, Real}} = nothing)

	"Plots multiple solution branches from a directory, using the closest available value for constant_param and automatically detecting which parameter varies. Optionally filters by varied_param_range = (min, max)."

	# validate directory
	if !isdir(dir_path)
		error("directory not found: $(dir_path)")
	end

	# collect all .jld2 files recursively
	files = String[]
	for (root, _, fs) in walkdir(dir_path)
		for f in fs
			endswith(f, ".jld2") && push!(files, joinpath(root, f))
		end
	end

	if isempty(files)
		error("no .jld2 files found in $(dir_path)")
	end

	# load all files and extract parameter values
	all_branches = Tuple{Matrix{Float64}, Constants}[]

	for file_path in files
		try
			solutions, constants, metadata = readSolution(file_path)
			push!(all_branches, (solutions, constants))
		catch
			# skip files that can't be read or don't have required fields
			continue
		end
	end

	if isempty(all_branches)
		error("no valid result files found in $(dir_path)")
	end

	# find the closest value of constant_param to the user's requested value
	param_values = [getfield(branch[2], constant_param) for branch in all_branches]
	unique_values = sort(unique(param_values))
	
	# find closest value
	closest_value = unique_values[argmin(abs.(unique_values .- constant_value))]
	
	# filter branches with this closest value
	matching_branches = Tuple{Matrix{Float64}, Constants}[]
	for (solutions, constants) in all_branches
		if getfield(constants, constant_param) == closest_value
			push!(matching_branches, (solutions, constants))
		end
	end

	# print what value is actually constant
	println("Plotting branches with $(constant_param) = $(closest_value)")

	# detect which parameter varies across the matching files
	# get all field names from the first constants struct (excluding z, dz)
	first_constants = matching_branches[1][2]
	param_names = filter(f -> f ∉ [:z, :dz, :N], fieldnames(typeof(first_constants)))
	
	# find parameters that vary
	varying_params = Symbol[]
	for param in param_names
		param == constant_param && continue  # skip the constant parameter
		
		values = [getfield(branch[2], param) for branch in matching_branches]
		unique_vals = unique(values)
		
		# if more than one unique value, this parameter varies
		if length(unique_vals) > 1
			push!(varying_params, param)
		end
	end

	if isempty(varying_params)
		@warn "no varying parameters detected - all branches have identical parameters"
		varying_param = :none
	elseif length(varying_params) == 1
		varying_param = varying_params[1]
	else
		# multiple parameters vary - use the first one and warn
		varying_param = varying_params[1]
		@warn "multiple parameters vary: $(varying_params), plotting by $(varying_param)"
	end

	# organize branches by the varying parameter
	branches_dict = Dict{Float64, Tuple{Matrix{Float64}, Constants}}()
	
	for (solutions, constants) in matching_branches
		if varying_param == :none
			# all branches identical, just plot them all
			key_val = rand()  # random key since they're all the same
		else
			key_val = getfield(constants, varying_param)
		end
		branches_dict[key_val] = (solutions, constants)
	end

	# sort by varying parameter value for consistent plotting
	sorted_values = sort(collect(keys(branches_dict)))

	# filter by varied_param_range if specified
	if !isnothing(varied_param_range) && varying_param != :none
		min_val, max_val = varied_param_range
		sorted_values = filter(v -> min_val <= v <= max_val, sorted_values)
		
		if isempty(sorted_values)
			error("no branches found within varied_param_range = $(varied_param_range)")
		end
		
		println("Filtered to $(length(sorted_values)) branches with $(varying_param) ∈ [$(min_val), $(max_val)]")
	end

	# create plot
	p = plot(legend = :best, size = (500, 500))
	
	for vary_val in sorted_values
		solutions, _ = branches_dict[vary_val]
		speeds = solutions[:, 1]
		a₁ = solutions[:, 3]
		
		# plot as line with label
		if varying_param == :none
			plot!(speeds, a₁, lw = 2, label = "branch")
		else
			plot!(speeds, a₁, 
				lw = 2, 
				label = "$(varying_param) = $(round(vary_val, digits=3))")
		end
	end

	xlabel!(L"c")
	ylabel!(L"a_1")

	return p
end

function plotBranchVariationsSubplots(dir_path::String; 
                                      constant_param::Symbol,
                                      constant_value::Real,
                                      varied_param_range::Union{Nothing, Tuple{Real, Real}} = nothing)

	"Creates a grid of subplots with one branch per subplot, where constant_param is held at the closest value to constant_value. Optionally filters by varied_param_range = (min, max)."

	# validate directory
	if !isdir(dir_path)
		error("directory not found: $(dir_path)")
	end

	# collect all .jld2 files recursively
	files = String[]
	for (root, _, fs) in walkdir(dir_path)
		for f in fs
			endswith(f, ".jld2") && push!(files, joinpath(root, f))
		end
	end

	if isempty(files)
		error("no .jld2 files found in $(dir_path)")
	end

	# load all files
	all_branches = Tuple{Matrix{Float64}, Constants}[]
	for file_path in files
		try
			solutions, constants, metadata = readSolution(file_path)
			push!(all_branches, (solutions, constants))
		catch
			continue
		end
	end

	if isempty(all_branches)
		error("no valid result files found in $(dir_path)")
	end

	# find the closest value of constant_param to the user's requested value
	param_values = [getfield(branch[2], constant_param) for branch in all_branches]
	unique_values = sort(unique(param_values))
	closest_value = unique_values[argmin(abs.(unique_values .- constant_value))]
	
	# filter branches with this closest value
	matching_branches = [branch for branch in all_branches 
	                     if getfield(branch[2], constant_param) == closest_value]
	
	println("Plotting branches with $(constant_param) = $(closest_value)")

	# detect which parameter varies
	first_constants = matching_branches[1][2]
	param_names = filter(f -> f ∉ [:z, :dz, :N], fieldnames(typeof(first_constants)))
	
	varying_params = Symbol[]
	for param in param_names
		param == constant_param && continue
		values = [getfield(branch[2], param) for branch in matching_branches]
		if length(unique(values)) > 1
			push!(varying_params, param)
		end
	end
	
	varying_param = isempty(varying_params) ? :none : varying_params[1]
	
	# organize branches by varying parameter
	branches_dict = Dict{Float64, Tuple{Matrix{Float64}, Constants}}()
	for (solutions, constants) in matching_branches
		key_val = varying_param == :none ? rand() : getfield(constants, varying_param)
		branches_dict[key_val] = (solutions, constants)
	end
	
	sorted_values = sort(collect(keys(branches_dict)))
	
	# filter by varied_param_range if specified
	if !isnothing(varied_param_range) && varying_param != :none
		min_val, max_val = varied_param_range
		sorted_values = filter(v -> min_val <= v <= max_val, sorted_values)
		
		if isempty(sorted_values)
			error("no branches found within varied_param_range = $(varied_param_range)")
		end
		
		println("Filtered to $(length(sorted_values)) branches with $(varying_param) ∈ [$(min_val), $(max_val)]")
	end
	
	# create one subplot per branch
	subplots = []
	for vary_val in sorted_values
		solutions, _ = branches_dict[vary_val]
		speeds = solutions[:, 1]
		a₁ = solutions[:, 3]
		
		# create subplot with title showing varying parameter value
		if varying_param == :none
			p = plot(speeds, a₁, lw = 4, legend = false, title = "branch",
			        xformatter = x -> string(round(x, digits=4)))
		else
			p = plot(speeds, a₁, lw = 4, legend = false, 
			        title = "$(varying_param) = $(round(vary_val, digits=3))",
			        xformatter = x -> string(round(x, digits=5)))
		end
		
		xlabel!(L"c")
		ylabel!(L"a_1")
		
		push!(subplots, p)
	end

	# determine grid layout (prefer roughly square)
	n_plots = length(subplots)
	n_cols = Int(ceil(sqrt(n_plots)))
	n_rows = Int(ceil(n_plots / n_cols))
	
	# create final plot
	p_final = plot(subplots..., layout = (n_rows, n_cols), size = (400 * n_cols, 400 * n_rows))
	
	return p_final
end

## LEGACY 

function plotting(solutions, index::Int, constants::Constants, shift_profiles = true)

	ii = index 

	# un-pack constants 
	z = constants.z
	L = constants.L

	branchN = length(solutions[:,1])

	# seperate coeffs and speeds
	coeffs = solutions[:,2:end]
	speeds = solutions[:,1];

	# create array for profiles
	profiles = zeros(branchN,length(z))
	
	# convert profiles
	Threads.@threads for i = 1:branchN
		profiles[i,:] .= fourierSeries(coeffs[i,:], z, L)[1]
	end

	# shift to the right profiles over by L
	if shift_profiles == true
		profiles = [profiles[:,Int(end/2)+1:end] profiles[:,1:Int(end/2)]]; nothing
	end

	# plot profiles 
	profile_plot = plot(z, profiles[ii,:], legend=false, title = "a1 = $(round(coeffs[ii,2], digits=3))", lw=2)
	xlabel!(L"z"); ylabel!(L"S")

	# plot coeffs 
	first_coeff = 0
	coeff_plot = scatter(abs.(coeffs[ii,first_coeff+1:end]), legend=false, title="k = $(round(ii*π/L))", xticks = :all, yaxis=:log)
	xlabel!("a$(first_coeff) to a$(length(coeffs[1,:])-1)")

	# plot branch
	branch_plot = scatter(speeds[1:ii], coeffs[1:ii,2], legend = false, markersize=4)
	xlabel!(L"c"); ylabel!(L"a_1")

	return profile_plot, branch_plot, coeff_plot
	
end

function plotting(solution_file::String)

	# load solutions
	solutions = readdlm("results/$(solution_file)/solutions_$(solution_file).dat")

	# load metadata
	meta = JSON.parsefile("results/$(solution_file)/meta_$(solution_file).json")

	# create constants
	constants = fuSimpleConstants(meta["N"], meta["L"], meta["b"], meta["λ2"], meta["vf"])

	profile_plot, branch_plot, coeff_plot = plotting(solutions, meta["branchN"], constants)
	convergence_plot = plot((meta["a1Vals"])[1:end-1], meta["iterations"], xlabel=L"a_1", ylabel="Iterations", legend=false, seriestype = :line, marker = :dot, markersize = 2)
	error_plot = plot((meta["a1Vals"])[1:end-1], meta["errors"], xlabel=L"a_1", ylabel="Error", title = "tol = $(meta["tol"])", legend=false, seriestype = :line, marker = :dot, markersize = 2, yaxis=:log10)
	# flag_plot = plot((meta["a1Vals"])[1:end-1], meta["flags"], xlabel=L"a_1", ylabel="Flag", legend=false, seriestype = :line, marker = :dot, markersize = 2)

	return profile_plot, branch_plot, coeff_plot, convergence_plot, error_plot
end