using Plots, LaTeXStrings

using SpecialFunctions
using Trapz
using LinearAlgebra
using Statistics

using DelimitedFiles, JSON

## Constants struct

struct Constants
	N::Int64  					# number of modes for solution S(z)
	L::Number                   # half-domain length
	
	# domain definition
	dz::Float64 				# domain spacing
    z::Vector{Float64} 			# domain vector of values (2N + 2 points)

	# physical geometry constants 
	b::Float64 					# inner rod radius
	λ2::Float64					# principal stretches

	vf::Float64					# fluid velocity 
	
	
	function Constants(N::Int64, L::Number, b::Float64, λ2::Float64, vf::Float64)
        dz = 2*L / (2*N+1)
        z = collect(-L:dz:L)

        new(N, L, dz, z, b, λ2, vf)
    end
end


## Main functions

# PROBLEM FUNCTIONS


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
		f(u::Vector{Float64}) = equations(u, constants, a1Vals[i], 1.0) 

		# solve for the current branch point + capture
		solutions[i,:], iterations[i], flags[i] = mySolver(f, initial_guess[i,:], tol = tol, solver = solver, max_iter = max_iter)

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
	# if directory exists, delete it
	if isdir("results/$(base_name)")
		rm("results/$(base_name)"; recursive = true)
	end

	# create directory to save results
	mkdir("results/$(base_name)")

	# save solutions, iterations and a1_vals to seperate files 
	writedlm("results/$(base_name)/solutions_$base_name.dat", solutions)

	# create a dictionary to save the solution metadata (constants + solver/solution parameters)
	meta = Dict( # Create a dictionary with the data
		"N" => constants.N,
		"L" => constants.L,
		"b" => constants.b,
		# "λ1" => constants.λ1,
		"λ2" => constants.λ2,
		"vf" => constants.vf,
		"tol" => tol,
		"branchN" => branchN,
		"solver" => solver,
		"max_iter" => max_iter,
		"initial_guess" => initial_guess[1,:],
		"a1Vals" => a1Vals,
		"iterations" => iterations,
		"errors" => errors, 
		"flags" => flags
	)

	# Write the dictionary to a JSON file
	open("results/$(base_name)/meta_$base_name.json", "w") do f
		JSON.print(f, meta)
	end

	return solutions, iterations

end

function c0(k, constants::Constants)
	## linearized wave speed for small amplitude waves c(k)

	b = constants.b
	λ2 = constants.λ2
	vf = constants.vf

	β0 = besseli(1, k) * besselk(1, k*b) - besseli(1, k*b) * besselk(1, k)
	β1 = besseli(1, k*b) * besselk(0, k) + besseli(0, k) * besselk(1, k*b) - (1/k)*β0

	# quadratic coeffs (Ac^2 + Bc + D) 
	A = k*β1 - 2*λ2^2 * β0 + β0
	B = 4 * vf * λ2 * β0
	D = -2 * vf^2 * β0

	# solve quadratic equation 
	c0 = solve_quadratic(A, B, D)

	return real(c0)
	
end

function β(n, k, b, S0)

	# i think this is wrong 

	beta1 = besseli.(1, complex.(k*b)) .* besselk.(n, complex.(k.*S0))
	beta2 = (-1)^n .* besselk.(1, complex.(k*b)) .* besseli.(n, complex.(k.*S0))

	return real.(beta1 .+ beta2)
end

function fourierSeries(coefficients::Vector{Float64}, domain, L::Number)
	
    N = length(coefficients) - 1
	
    S = zeros(length(domain))  		# profile S
    Sz = zeros(length(domain))  	# first derivative Sz	
    Szz = zeros(length(domain))  	# second derivative Szz

    Threads.@threads for i in 1:Int(length(domain))

		x = domain[i]
		
        # Calculate the series and its derivatives at each point x
        S[i] = coefficients[1] 
		
        for n in 1:N

			k = n*π/L

            S[i] += coefficients[n + 1] * cos(k * x)
            Sz[i] -= k * coefficients[n + 1] * sin(k * x)
            Szz[i] -= k^2 * coefficients[n + 1] * cos(k * x)
        end
		
    end

    return S, Sz, Szz
end

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

function recompute_solutions(file_path::String)

	# get names of the files in the results folder
	files = readdir(file_path)

	# loop over the files
	for file in files
		
		# get the parameters from the file name 
		# (note that the second parameter is the tolerance, which has a decimal point in it, so we need to add the "1" to the start of the string to parse it as a float)
		parameters = split(file, ".")
		N = parse(Int, parameters[1])
		tol = parse(Float64, "1." * parameters[3])
		branchN = parse(Int, parameters[4])
		solver = parameters[5]

		# print each for the first 10 files 
		if file == files[1:10]
			println("N: ", N)
			println("tol: ", tol)
			println("branchN: ", branchN)
			println("solver: ", solver)
		end

		# initialize everything
		a1Vals = collect(range(0.001, 0.33, branchN + 1))		  # a1 values

		constants = Constants(N,π,1.5,0.1) 					      # constants  
		k1 = 1*π/constants.L			   					      # wave number
		cInitial = c0(k1, constants);	   					      # wave speed 

		initial_guess = (1e-16).*ones(branchN+1, constants.N + 2) # initial guess 
		initial_guess[1,1:4] = [cInitial, 1.0, a1Vals[1], 1e-16]

		# solve 
		bifurcation(initial_guess, a1Vals, branchN, constants, tol = tol, solver = Symbol(solver), max_iter = 40000)

	end

end

function solve_quadratic(a, b, c)
    Δ = b^2 - 4*a*c
    if Δ < 0
		print("Solution (c0?) is complex!")
        return Complex.((-b + sqrt(Δ)) / (2*a))  # Return complex roots
    else
        return (-b + sqrt(Δ)) / (2*a)  # Return real roots
    end
end