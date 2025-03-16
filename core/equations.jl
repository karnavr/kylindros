"General equations for the cylindrical AFM formulation for an arbitrary wall model."

function equations(unknowns, constants::Constants, a₁::Float64, a₀::Float64)

	# Returns the N + 2 equations that we want to solve for:

	# problem constants 
	z = constants.z::Vector{Float64}
	N = constants.N::Int64
	b = constants.b::Float64
	L = constants.L::Number
	λ2 = constants.λ2::Float64
	vf = constants.vf::Float64

	c = unknowns[1]
	coeffs = unknowns[2:N+2] # N + 1 coeffs

	a0 = coeffs[1]
	a1 = coeffs[2]

	S, Sz, _ = fourierSeries(coeffs, z, L)

	# Use the element type of unknowns for the output arrays
	T = eltype(unknowns)
	integrands = zeros(T, N, length(z))  # Use the same type as unknowns
	integrals = zeros(T, N)              # Use the same type as unknowns
	eqs = zeros(T, N+2)                  # Use the same type as unknowns

	# define common factor in equations 
	Szsq = 1 .+ (Sz.^2);

	# define wall model
	w = wall_model(constants, c, S)
	one_p = Szsq .* (c^2 .- 2 .* w)

	Threads.@threads for n = 1:N

		k = n*π/L 
		
	    one = k .* S .* sqrt.(Complex.(one_p))
	    two = besselk.(1, Complex.(k * b)) .* besseli.(1, Complex.(k .* S)) .- besseli.(1, Complex.(k * b)) .* besselk.(1, Complex.(k .* S))

	    integrands[n, :] = real.(one .* two)
		
	    # Normalize the integrand before integration to prevent numerical issues
	    integrands[n, :] ./= maximum(abs.(integrands[n, :]))
	    integrals[n] = trapz(z, integrands[n, :] .* cos.(k .* z))

	end

	eqs[1:N] = real.(integrals)
	eqs[N+1] = abs(a0 - a₀)
	eqs[N+2] = abs(a1 - a₁)

	return eqs
	
end
