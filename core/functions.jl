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