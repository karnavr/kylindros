"General equations for the cylindrical AFM formulation with an arbitrary wall model."

function equations(unknowns, constants::Constants, a₁::Float64, a₀::Float64)

	# Returns the N + 2 equations that we want to solve for:
	# - N integrals 
	# - 2 extra equations setting values of a0 and a1

	# problem constants (unpack and assign to variables)
	unpackConstants(constants)

	c = unknowns[1]
	coeffs = unknowns[2:N+2] # N + 1 coeffs

	a0 = coeffs[1]
	a1 = coeffs[2]

	S, Sz, Szz = fourierSeries(coeffs, z, L)

	# Use the element type of unknowns for the output arrays
	T = eltype(unknowns)
	integrands = zeros(T, N, length(z))
	integrals = zeros(T, N)
	eqs = zeros(T, N+2)            

	# define common factor in equations 
	Szsq = 1 .+ (Sz.^2);

	# define wall model
	if constants isa ferrofluidConstants
		one_p = (Szsq).*((c^2)./2 .- 1 ./ (S.*sqrt.(Szsq)) .+ Szz./(Szsq.^(3/2)) .+ B./(2 .* S.^2) .+ E);
	else 
		w = wall_model(constants, c, S)
		one_p = Szsq .* (c^2 .- 2 .* w)
	end
	

	Threads.@threads for n = 1:N

		k = n*π/L 
		
	    one = k .* S .* sqrt.(Complex.(one_p))
	    # two = besselk.(1, Complex.(k * b)) .* besseli.(1, Complex.(k .* S)) .- besseli.(1, Complex.(k * b)) .* besselk.(1, Complex.(k .* S))
		two = β_scaled(k, S, b)

	    integrands[n, :] = real.(one .* two)
		
	    # Normalize the integrand before integration to prevent numerical issues
	    integrands[n, :] ./= maximum(abs.(integrands[n, :]))
	    integrals[n] = trapz(z, integrands[n, :] .* cos.(k .* z))

	end

	# Threads.@threads for n = 1:N
	# 	k = n*π/L
	# 	y = real.(k .* S .* sqrt.(Complex.(one_p)) .* β_scaled(k, S, b))
	# 	y ./= maximum(abs.(y))
	# 	integrals[n] = trapz(z, y .* cos.(k .* z))
	# end

	eqs[1:N] = real.(integrals)
	eqs[N+1] = abs(a0 - a₀)
	eqs[N+2] = abs(a1 - a₁)

	return eqs
	
end

function β_scaled(k, S, b)
    
    kb = k*b
    kS = k.*S

    # preallocate L₁ and L₂
    L₁ = zeros(length(S))
    L₂ = zeros(length(S))

    # compute L₁ and L₂
    L₁ .= kS .- kb .+ log.(SpecialFunctions.besselkx.(1, Complex.(kb))) .+ log.(SpecialFunctions.besselix.(1, Complex.(kS)))
    L₂ .= kb .- kS .+ log.(SpecialFunctions.besselix.(1, Complex.(kb))) .+ log.(SpecialFunctions.besselkx.(1, Complex.(kS)))

    # preallocate β
    β = zeros(length(S))

    # compute max of L₁ and L₂ (ensure stabilization is applied for each element in domain)
    M = max.(L₁, L₂)

    # compute β
    β .= exp.(M) .* (exp.(L₁ .- M) .- exp.(L₂ .- M))

    return real.(β)
end

# function β_scaled(k, S, b)
    
#     kb = k*b
#     kS = k.*S

#     # preallocate L₁ and L₂
#     L₁ = zeros(ComplexF64, length(S))
#     L₂ = zeros(ComplexF64, length(S))

#     # Safely compute Bessel functions with explicit complex conversion
#     # and handle potential numerical issues
#     for i in eachindex(S)
#         # Ensure positive arguments when possible to avoid unnecessary complex calculations
#         if kS[i] > 0 && kb > 0
#             # Real branch - no need for complex conversion
#             L₁[i] = kS[i] - kb + log(SpecialFunctions.besselkx(1, kb)) + log(SpecialFunctions.besselix(1, kS[i]))
#             L₂[i] = kb - kS[i] + log(SpecialFunctions.besselix(1, kb)) + log(SpecialFunctions.besselkx(1, kS[i]))
#         else
#             # Handle potentially complex branch carefully
#             L₁[i] = kS[i] - kb + log(SpecialFunctions.besselkx(1, Complex(abs(kb)))) + log(SpecialFunctions.besselix(1, Complex(abs(kS[i]))))
#             L₂[i] = kb - kS[i] + log(SpecialFunctions.besselix(1, Complex(abs(kb)))) + log(SpecialFunctions.besselkx(1, Complex(abs(kS[i]))))
#         end
#     end

#     # preallocate β
#     β = zeros(ComplexF64, length(S))

#     # compute max of real parts of L₁ and L₂ for numerical stability
#     M = max.(real.(L₁), real.(L₂))

#     # compute β with numerical stabilization
#     for i in eachindex(S)
#         β[i] = exp(M[i]) * (exp(L₁[i] - M[i]) - exp(L₂[i] - M[i]))
        
#         # Safety check - if imaginary part is significant, something is wrong
#         if abs(imag(β[i])) > 1e-10 * abs(real(β[i]))
#             @warn "Large imaginary component detected at index $i: $(β[i])"
#         end
#     end

#     # Return only the real part, which should be the physically meaningful result
#     return real.(β)
# end

# function β_scaled(k, S, b)
#     kb = k*b
#     kS = k.*S
    
#     # Define threshold for "large" arguments where direct computation is unsafe
#     LARGE_ARG_THRESHOLD = 700.0  # Typical threshold for numerical stability
    
#     # Preallocate arrays
#     L₁ = zeros(ComplexF64, length(S))
#     L₂ = zeros(ComplexF64, length(S))
#     β = zeros(ComplexF64, length(S))
    
#     for i in eachindex(S)
#         # Handle each point carefully with bounds checking
#         try
#             # Compute L₁ and L₂ safely
#             if abs(kb) > LARGE_ARG_THRESHOLD || abs(kS[i]) > LARGE_ARG_THRESHOLD
#                 # Use asymptotic approximations for large arguments
#                 # For large x: besselkx(1,x) * exp(x) ≈ sqrt(π/(2x))
#                 # For large x: besselix(1,x) * exp(-x) ≈ 1/sqrt(2πx)
                
#                 # logarithm of approximations for numerical stability
#                 if abs(kb) > LARGE_ARG_THRESHOLD
#                     log_besselkx_kb = 0.5*log(π/(2*abs(kb)))
#                     log_besselix_kb = -0.5*log(2*π*abs(kb))
#                 else
#                     # For moderate values, use actual functions
#                     log_besselkx_kb = log(SpecialFunctions.besselkx(1, Complex(kb)))
#                     log_besselix_kb = log(SpecialFunctions.besselix(1, Complex(kb)))
#                 end
                
#                 if abs(kS[i]) > LARGE_ARG_THRESHOLD
#                     log_besselkx_kS = 0.5*log(π/(2*abs(kS[i])))
#                     log_besselix_kS = -0.5*log(2*π*abs(kS[i]))
#                 else
#                     # For moderate values, use actual functions
#                     log_besselkx_kS = log(SpecialFunctions.besselkx(1, Complex(kS[i])))
#                     log_besselix_kS = log(SpecialFunctions.besselix(1, Complex(kS[i])))
#                 end
                
#                 # Compute L₁ and L₂ using the (log of) approximations
#                 L₁[i] = kS[i] - kb + log_besselkx_kb + log_besselix_kS
#                 L₂[i] = kb - kS[i] + log_besselix_kb + log_besselkx_kS
#             else
#                 # Standard computation for moderate arguments
#                 L₁[i] = kS[i] - kb + log(SpecialFunctions.besselkx(1, Complex(kb))) + 
#                         log(SpecialFunctions.besselix(1, Complex(kS[i])))
#                 L₂[i] = kb - kS[i] + log(SpecialFunctions.besselix(1, Complex(kb))) + 
#                         log(SpecialFunctions.besselkx(1, Complex(kS[i])))
#             end
            
#             # Numerical stabilization - use max of real parts
#             M = max(real(L₁[i]), real(L₂[i]))
            
#             # Compute β with stabilization
#             β[i] = exp(M) * (exp(L₁[i] - M) - exp(L₂[i] - M))
            
#         catch e
#             # Handle any exceptions during computation
#             if isa(e, DomainError) || isa(e, AmosException)
#                 # If we get a domain error or Amos exception, use the asymptotic limit
#                 # For very large or small values, L₁ and L₂ will tend to cancel out
#                 β[i] = 0.0
#                 @warn "Exception in β_scaled calculation at index $i with kS=$(kS[i]), kb=$kb: $e. Setting result to 0."
#             else
#                 # Rethrow other errors
#                 rethrow(e)
#             end
#         end
#     end
    
#     # Return real part as the physically meaningful result
#     return real.(β)
# end