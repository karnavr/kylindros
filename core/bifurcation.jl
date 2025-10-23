"Solves for a bifurcation branch of solutions."

using LinearAlgebra

function bifurcation(initial_guess::Matrix{Float64}, a1Vals, branchN::Int64, constants::Constants; tol = 1e-12, solver = :myNewtonRaphson, max_iter = 1000, overwrite = false, save_dir::String = "", verbose = true)::Tuple{Matrix{Float64}, Constants, Dict{String, Any}}

	"Compute the bifurcation branch for branchN branch points and provided a₁ values, starting at the given initial guess, using a continuation scheme."

	# check type of the a1Vals argument 
	if typeof(a1Vals) != Vector{Float64}

		# create a vector of length branchN starting from 0.001 with stepsize a1Vals
		a1Vals = collect(range(0.001, step = a1Vals, length = branchN + 1))

	end


	# check for existing solution and determine save directory
	existing_result, dir, model_name = handleExistingSolution(constants, branchN, a1Vals, tol, solver, save_dir, overwrite)
	
	if !isnothing(existing_result)
		return existing_result
	end
	
	# initialize solution array
	solutions = zeros(Float64, branchN, constants.N+2)

	# initialize number of iterations array
	iterations = zeros(Int64, branchN)

	# initialize flags array (for each branch point)
	flags = fill(NaN, branchN)

	# condition numbers per branch point (only filled for :NLSolver for now)
	condition_numbers = zeros(Float64, branchN)

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

			# compute Jacobian condition number at converged solution
			funcJ(u) = equations(u, constants, a1Vals[i], 1.0)
			J = finite_diff_jacobian(funcJ, sol.u)
			condition_numbers[i] = cond(J)

		else

			# define the set of equations/function to solve: f(x) = 0
			func(u) = equations(u, constants, a1Vals[i], 1.0)

			# solve for the current branch point + capture
			solutions[i,:], iterations[i], flags[i] = mySolver(func, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)
			
		end

		# update intial guess with current solution
		initial_guess[i+1,:] = solutions[i,:]

        # print progress for every 10% of branch points 
        if i % Int(round(0.1*branchN)) == 0 && verbose
                println("Branch point $i of $branchN, $(Int(iterations[i])) iterations.")
        end

		# check the average values of the last 20% of coefficients every 5% of branch points (TODO: maybe there is a more precise way to do this)
		# and if it's above 1e-15, zero the last 20% of coefficients
		if i % max(1, Int(round(0.02*branchN))) == 0 # max() to ensure it's always greater than 1
			if mean(abs.(initial_guess[i, end - Int(round(0.2*length(initial_guess[1,:]))):end])) > 1e-15
				initial_guess[i+1, end - Int(round(0.2*length(initial_guess[1,:]))):end] .= 0
			end
		end

	end

	# end timer and compute time 
	computation_time = time() - start_time

	# compute error (L2 norm) for each branch point
	errors = zeros(Float64, branchN)
	for i = 1:branchN
		errors[i] = norm(equations(solutions[i,:], constants, a1Vals[i], 1.0))
	end


	## SAVING RESULTS

	metadata = Dict{String, Any}(
		"model" => model_name,
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
		"condition_numbers" => condition_numbers,
	)

	# build filename and create directory if needed
	filename = generateFileName(metadata, constants)
	mkpath(dir)

	# save the components directly at the top level
	@save joinpath(dir, "$(filename).jld2") solutions constants metadata
	if verbose
		println("Saved solution branch to $(filename).jld2")
	end

	return solutions, constants, metadata

end