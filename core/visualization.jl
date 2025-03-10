"Plotting functions for the results of the numerical continuation."

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
	constants = Constants(meta["N"], meta["L"], meta["b"], meta["λ2"], meta["vf"])

	profile_plot, branch_plot, coeff_plot = plotting(solutions, meta["branchN"], constants)
	convergence_plot = plot((meta["a1Vals"])[1:end-1], meta["iterations"], xlabel=L"a_1", ylabel="Iterations", legend=false, seriestype = :line, marker = :dot, markersize = 2)
	error_plot = plot((meta["a1Vals"])[1:end-1], meta["errors"], xlabel=L"a_1", ylabel="Error", title = "tol = $(meta["tol"])", legend=false, seriestype = :line, marker = :dot, markersize = 2, yaxis=:log10)
	# flag_plot = plot((meta["a1Vals"])[1:end-1], meta["flags"], xlabel=L"a_1", ylabel="Flag", legend=false, seriestype = :line, marker = :dot, markersize = 2)

	return profile_plot, branch_plot, coeff_plot, convergence_plot, error_plot
end
