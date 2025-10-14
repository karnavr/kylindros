"Solves for a bifurcation branch of solutions."

using LinearAlgebra

function bifurcation(initial_guess::AbstractMatrix{T}, a1Vals::AbstractVector{T}, branchN::Int64, constants::Constants; tol::T = T(1e-12), solver = :myNewtonRaphson, max_iter = 1000, overwrite = false, save_subdir::String = "", verbose = true) where {T}

	"Compute the bifurcation branch for branchN branch points and provided a₁ values, starting at the given intial guess"

	# check if the solution branch already exists
	existing_filename = solutionExists(constants, branchN, float(tol); subdir = save_subdir)

	# if the solution branch already exists, return the existing solution or delete it based on overwrite flag
	if existing_filename != false
		model_name = getModelName(constants)
		dir = joinpath(results_dir(save_subdir), model_name)
		full_existing_filename = joinpath(dir, string(existing_filename))
		
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
	solutions = zeros(T, branchN, constants.N+2)

	# initialize convergence array
	iterations = zeros(Int, branchN)

	# initialize flags array (for each branch point)
	flags = zeros(Int, branchN)

	# condition numbers per branch point (only filled for :NLSolver for now)
	condition_numbers = zeros(Float64, branchN)

	# start timer
	start_time = time()
	
	for i = 1:branchN

		
		if solver == :NLSolver

			# define the set of equations/function to solve: f(x) = 0
			f(u, p) = equations(u, constants, a1Vals[i], T(1))

			# define the problem
			problem = NonlinearProblem(f, initial_guess[i,:])

			# solve the problem with more conservative settings
			sol = solve(problem, NewtonRaphson(; autodiff = AutoFiniteDiff()))
			# sol = solve(problem, LevenbergMarquardt(; autodiff = AutoFiniteDiff()))

			solutions[i,:] = sol.u
			iterations[i] = sol.stats.nsteps
			# flags[i] = sol.stats.status

			# compute Jacobian condition number at converged solution
			funcJ(u) = equations(u, constants, a1Vals[i], T(1))
			J = finite_diff_jacobian(funcJ, sol.u)
			condition_numbers[i] = cond(J)

		else

			# define the set of equations/function to solve: f(x) = 0
			func(u) = equations(u, constants, a1Vals[i], T(1))

			# solve for the current branch point + capture
			solutions[i,:], iterations[i], flags[i] = mySolver(func, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)
			
		end

		# update intial guess with current solution
		initial_guess[i+1,:] = solutions[i,:]

        # print progress for every 10% of branch points 
        if i % Int(round(0.01*branchN)) == 0 && verbose
                println("Branch point $i of $branchN, $(Int(iterations[i])) iterations.")
        end

		# check the average values of the last 20% of coefficients every 5% of branch points (TODO: maybe there is a more precise way to do this)
		# and if it's above 1e-15, zero the last 20% of coefficients
		if i % max(1, Int(round(0.02*branchN))) == 0 # max() to ensure it's always greater than 1
			if mean(abs.(initial_guess[i, end - Int(round(0.2*length(initial_guess[1,:]))):end])) > 1e-15
				initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= T(0)
			end
		end

	end

	# end timer and compute time 
	computation_time = time() - start_time

	# compute error (L2 norm) for each branch point
	errors = zeros(T, branchN)
	Threads.@threads for i = 1:branchN
		errors[i] = norm(equations(solutions[i,:], constants, a1Vals[i], T(1)))
	end


	## SAVING RESULTS

	metadata = Dict(
		"model" => getModelName(constants),
		"timestamp" => getCurrentDateToString(),
		"computation_time" => computation_time,
		"solver" => solver,
		"max_iter" => max_iter,
		"tol" => float(tol),
		"branchN" => branchN,
		"a1Vals" => float.(a1Vals),
		"iterations" => iterations,
		"errors" => float.(errors),
		"flags" => flags,
		"condition_numbers" => condition_numbers,
	)

	# build filename and destination directory
	filename = string(generateFileName(metadata), "_N", constants.N)
	model_name = metadata["model"]
	dir = joinpath(results_dir(save_subdir), model_name)
	mkpath(dir)

	# Save the components directly at the top level
	@save joinpath(dir, string(filename, ".jld2")) solutions constants metadata
	if verbose
		println("Saved solution branch to $(filename).jld2")
	end

	return solutions, constants, metadata

end
