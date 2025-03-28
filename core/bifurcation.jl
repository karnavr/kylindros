"Solves for a bifurcation branch of solutions."

function bifurcation(initial_guess, a1Vals, branchN::Int64, constants::Constants; tol = 1e-12, solver = :myNewtonRaphson, max_iter = 1000, overwrite = false)

	"Compute the bifurcation branch for branchN branch points and provided a₁ values, starting at the given intial guess"

	# check type of the a1Vals argument 
	if typeof(a1Vals) != Vector{Float64}

		# create a vector of length branchN starting from 0.001 with stepsize a1Vals
		a1Vals = collect(range(0.001, step = a1Vals, length = branchN + 1))

	end


	# check if the solution branch already exists
	existing_filename = solutionExists(constants, branchN, tol)

	# if the solution branch already exists, return the existing solution or delete it based on overwrite flag
	if existing_filename != false
		model_name = getModelName(constants)
		full_existing_filename = "results/$(model_name)/$(existing_filename)"
		
		if !overwrite
			println("Solution branch already exists.")
			println("File: $(full_existing_filename)")
			return load(full_existing_filename)["solutions"], load(full_existing_filename)["constants"], load(full_existing_filename)["metadata"]
		else
			println("Solution branch already exists, but overwrite flag is set.")
			println("Deleting existing file: $(full_existing_filename)")
			rm(full_existing_filename)
		end
	end
	
	# initialize solution array
	solutions = zeros(branchN, constants.N+2)

	# initialize convergence array
	iterations = zeros(branchN)

	# initialize flags array
	flags = randn(branchN)

	# start timer
	start_time = time()
	
	for i = 1:branchN

		
		if solver == :NLSolver

			# define the set of equations/function to solve: f(x) = 0
			f(u, p) = equations(u, constants, a1Vals[i], 1.0)

			# define the problem
			problem = NonlinearProblem(f, initial_guess[i,:])

			# solve the problem with more conservative settings
			sol = solve(problem, NewtonRaphson(; autodiff = AutoFiniteDiff()))
			# sol = solve(problem, LevenbergMarquardt(; autodiff = AutoFiniteDiff()))

			solutions[i,:] = sol.u
			iterations[i] = sol.stats.nsteps
			# flags[i] = sol.stats.status

		else

			# define the set of equations/function to solve: f(x) = 0
			func(u) = equations(u, constants::Constants, a1Vals[i], 1.0)

			# solve for the current branch point + capture
			solutions[i,:], iterations[i], flags[i] = mySolver(func, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)
			
		end

		# update intial guess 
		initial_guess[i+1,:] = solutions[i,:]

        # print progress for every 10% of branch points 
        if i % Int(round(0.1*branchN)) == 0
            println("Branch point $i of $branchN, $(Int(iterations[i])) iterations.")
        end

		# check the average values of the last 20% of coefficients every 5% of branch points
		# and if it's above 1e-15, zero the last 20% of coefficients
		if i % max(1, Int(round(0.02*branchN))) == 0 # max() to ensure it's always greater than 1
			if mean(abs.(initial_guess[i, end - Int(round(0.2*length(initial_guess[1,:]))):end])) > 1e-15
				initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= 0
			end
		end

	end

	# end timer and compute time 
	computation_time = time() - start_time

	# compute error for each branch point
	errors = zeros(branchN)

	Threads.@threads for i = 1:branchN
		errors[i] = norm(equations(solutions[i,:], constants, a1Vals[i], 1.0))
	end


	## SAVING RESULTS

	metadata = Dict(
		"model" => getModelName(constants),
		"timestamp" => getCurrentDateToString(),
		"computation_time" => computation_time,
		"solver" => solver,
		"max_iter" => max_iter,
		"tol" => tol,
		"branchN" => branchN,
		"a1Vals" => a1Vals,
		"iterations" => iterations,
		"errors" => errors,
		"flags" => flags,
	)

	filename = generateFileName(metadata)

	# Save the components directly at the top level
	@save "../results/$(metadata["model"])/$(filename).jld2" solutions constants metadata
	println("Saved solution branch to $(filename).jld2")

	return solutions, constants, metadata

end