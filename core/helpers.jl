function β(n, k, b, S0)

	# i think this is wrong 

	beta1 = besseli.(1, complex.(k*b)) .* besselk.(n, complex.(k.*S0))
	beta2 = (-1)^n .* besselk.(1, complex.(k*b)) .* besseli.(n, complex.(k.*S0))

	return real.(beta1 .+ beta2)
end

function fourierSeries(coefficients::Vector{Float64}, domain::Vector{Float64}, L::Real)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
	
    N = length(coefficients) - 1
	n_domain = length(domain)
	
    S = zeros(Float64, n_domain)  	# profile S
    Sz = zeros(Float64, n_domain)  	# first derivative Sz	
    Szz = zeros(Float64, n_domain)  	# second derivative Szz

	k_values = [(n*π/L) for n in 1:N]
	
	# Add constant offset to S
	S .= coefficients[1] 
	
	for n in 1:N
		k = k_values[n]
		coeff = coefficients[n + 1]

		@. S += coeff * cos(k * domain)
		@. Sz -= k * coeff * sin(k * domain)
		@. Szz -= k^2 * coeff * cos(k * domain)
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

function solve_quadratic(a::AbstractArray{Float64}, b::AbstractArray{Float64}, c::AbstractArray{Float64})
    Δ = b.^2 - 4 .* a .* c
    if any(Δ .< 0)
	print("Solution (c0?) is complex!")
        return Complex.((-b .+ sqrt.(Complex.(Δ))) ./ (2 .* a))  # Return complex roots
    else
        return (-b .+ sqrt.(Complex.(Δ))) ./ (2 .* a)  # Return real roots
    end
end

function solve_quadratic(a::Float64, b::Float64, c::Float64)
    Δ = b^2 - 4*a*c
    if Δ < 0
	print("Solution (c0?) is complex!")
        return Complex.((-b + sqrt(Δ)) / (2*a))  # Return complex roots
    else
        return (-b + sqrt(Δ)) / (2*a)  # Return real roots
    end
end