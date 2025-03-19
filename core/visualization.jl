"Plotting functions for the results of the numerical continuation."

## Functions for each type of plot

function plot_profiles(solutions, constants::Constants; shift_profiles = true)

	## plot three profiles from the solutions branch (equally spaced in the last 0.75 of the branch)

	# get needed constants
	L = constants.L
	z = constants.z
	branchN = length(solutions[:,1])

	# seperate coeffs and speeds
	coeffs = solutions[:,2:end]

	# pick three equally spaced indices from the last 0.75 of the branch
	indices = [Int(round(0.5*branchN)), Int(round(0.75*branchN)), branchN]

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
		plot!(z, profiles[index,:], label = "a₁ = $(round(solutions[index,3], digits=3))", lw=2, linestyle = linestyles[i])
	end
	xlabel!(L"z"); ylabel!(L"S")

	return p
end

function plot_branch(solutions, metadata)

	## plot the branch of solutions

	# get errors
	errors = metadata["errors"]

	# seperate coeffs and speeds
	speeds = solutions[:,1]
	a₁ = solutions[:,3]

	# plot the branch with error-based coloring
	p = scatter(speeds, a₁, 
		marker_z = log10.(abs.(errors)),         # color by log10 of errors
		color = :plasma,           # colormap choice
		colorbar_title = "log₁₀(Error)",  # colorbar label
		colorbar_titlefontsize = 6,
		colorbar_tickfontsize = 6,
		colorbar = false,
		legend = false, 
		size = (750,500))
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
		plot!(abs.(coeffs[index,:]), 
			yaxis=:log, 
			label = "a₁ = $(round(solutions[index,3], digits=3))", 
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

	branchN = metadata["branchN"]

	profileplot = plot_profiles(solutions, constants, shift_profiles = true)
	branchplot = plot_branch(solutions, metadata)
	coeffplot = plot_coeffs(solutions, 1:(branchN/10):branchN)
	errorplot = plot_error(solutions, metadata)

	p = plot(profileplot, branchplot, coeffplot, errorplot, layout = (2,2), size = (1000,1000))

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
		speeds[i,:] .= c0(k_range, constants)
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