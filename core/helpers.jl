function β(n, k, b, S0)

	# i think this is wrong 

	beta1 = besseli.(1, complex.(k*b)) .* besselk.(n, complex.(k.*S0))
	beta2 = (-1)^n .* besselk.(1, complex.(k*b)) .* besseli.(n, complex.(k.*S0))

	return real.(beta1 .+ beta2)
end

function fourierSeries(coefficients::Vector{Float64}, domain::Vector{Float64}, L::Real)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
	
    N = length(coefficients) - 1
	n_domain = length(domain)
	
    S = zeros(Float64, n_domain)
    Sz = zeros(Float64, n_domain)
    Szz = zeros(Float64, n_domain)
    Szzz = zeros(Float64, n_domain)
    Szzzz = zeros(Float64, n_domain)

	# Add constant offset to S
	S .= coefficients[1] 
	
	# Pre-allocate buffers for trig values to avoid recomputation
	cos_kz = Vector{Float64}(undef, n_domain)
	sin_kz = Vector{Float64}(undef, n_domain)
	
	π_over_L = π / Float64(L)  # Compute once and ensure Float64
	
	for n in 1:N
		k = n * π_over_L  # Compute inline, no array allocation
		coeff = coefficients[n + 1]
		
		# Compute trig values once, reuse for all derivatives
		@. cos_kz = cos(k * domain)
		@. sin_kz = sin(k * domain)
		
		@. S += coeff * cos_kz
		@. Sz -= k * coeff * sin_kz
		@. Szz -= (k^2) * coeff * cos_kz
		@. Szzz += (k^3) * coeff * sin_kz
		@. Szzzz += (k^4) * coeff * cos_kz
	end
		
    return S, Sz, Szz, Szzz, Szzzz
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

function compute_cauchy_error(fileN_path::String, fileNm2_path::String, idx::Int)

	"Compute Cauchy convergence error E_N = ||S_N - S_{N-2}||_{L2} for branch index idx. Returns (E_N, a1)."

	# ========== Load solutions ==========
	solsN, constsN, _ = readSolution(fileN_path)
	solsNm2, constsNm2, _ = readSolution(fileNm2_path)

	# ========== Validation checks ==========
	if constsN.L != constsNm2.L
		error("L mismatch between N=$(constsN.N) and N=$(constsNm2.N)")
	end

	if idx > size(solsN, 1) || idx > size(solsNm2, 1)
		error("Index $idx out of bounds (N branch: $(size(solsN,1)), N-2 branch: $(size(solsNm2,1)))")
	end

	# ========== Extract Fourier coefficients ==========
	# Extract from column 2 onward (includes amplitude A and Fourier coefficients)
	coeffsN = solsN[idx, 2:end]
	coeffsNm2 = solsNm2[idx, 2:end]

	# ========== Reconstruct profiles on common grid ==========
	# Use finer grid from N run for evaluation
	z_grid = constsN.z
	L = constsN.L

	S_N = fourierSeries(coeffsN, z_grid, L)[1]
	S_Nm2 = fourierSeries(coeffsNm2, z_grid, L)[1]

	# ========== Compute continuous L2 norm ==========
	diff_sq = (S_N .- S_Nm2).^2
	E_N = sqrt(trapz(z_grid, diff_sq))

	# Extract a1 coefficient for labeling (column 3 of solution array)
	a1 = solsN[idx, 3]

	return E_N, a1
end