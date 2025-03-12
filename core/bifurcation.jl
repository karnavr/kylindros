"Solves for a bifurcation branch of solutions."

function bifurcation(initial_guess::Matrix{Float64}, a1Vals, branchN::Int64, constants::Constants; tol = 1e-8, solver = :NewtonRaphson, max_iter = 1000)

	## compute the bifurcation branch for branchN branch points and provided a₁ values, starting at the given intial guess

	# check type of the a1Vals argument 
	if typeof(a1Vals) != Vector{Float64}

		# create a vector of length branchN starting from 0.001 with stepsize a1Vals
		a1Vals = collect(range(0.001, step = a1Vals, length = branchN + 1))

	end

	# create base file name (num_modes).(tolerance).(branchN).(solver).(max_iter)
	base_name = "$(constants.N).$(tol).$(branchN).$(solver)"
	println("$(base_name)")

	# check if the directory exists -> if it does, display a warning
	if isdir("results/$(base_name)")
		println("Warning: Solution branch already exists. Overwriting files.")
		# NOTE: will have to eventually add a check for all the constants, etc.
	end
	
	# initialize solution array
	solutions = zeros(branchN, constants.N+2)

	# initialize convergence array
	iterations = zeros(branchN)

	# initialize flags array
	flags = randn(branchN)
	
	for i = 1:branchN

		# define the set of equations/function to solve: f(x) = 0
		f(u::Vector{Float64}) = equations(u::Vector{Float64}, constants::Constants, a1Vals[i], 1.0) 

		# solve for the current branch point + capture
		solutions[i,:], iterations[i], flags[i] = mySolver(f::Function, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)

		# update intial guess 
		initial_guess[i+1,:] = solutions[i,:]

        # print progress for every 10% of branch points 
        if i % Int(round(0.1*branchN)) == 0
            println("Branch point $i of $branchN, $(Int(iterations[i])) iterations.")
			# initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= 0
        end

		# check the average values of the last 20% of coefficients every 1% of branch points
		# and if it's above 1e-15, zero the last 20% of coefficients
		if i % max(1, Int(round(0.01*branchN))) == 0 # max() to ensure it's always greater than 1
			if mean(abs.(initial_guess[i, end - Int(round(0.2*length(initial_guess[1,:]))):end])) > 1e-15
				initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= 0
			end
		end

	end

	# compute error for each branch point
	errors = zeros(branchN)

	Threads.@threads for i = 1:branchN
		errors[i] = norm(equations(solutions[i,:], constants, a1Vals[i], 1.0))
	end


	## SAVING RESULTS

	metadata = Dict(
		"model" => getModelName(constants),
		"solver" => solver,
		"max_iter" => max_iter,
		"tol" => tol,
		"branchN" => branchN,
		"a1Vals" => a1Vals,
		"iterations" => iterations,
		"errors" => errors,
		"flags" => flags
	)

	results_dict = createResultsDict(constants, solutions, metadata)
	filename = generateFileName(results_dict)

	@save "results/$(filename).jld2" results_dict


	return solutions, iterations

end
