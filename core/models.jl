"Constants structs, linearized wave speed and w term for each wall model."

abstract type Constants end

## FERROFLUID (BLYTH & PARAU)

struct ferrofluidConstants <: Constants
	N::Int64  					# number of modes for solution S(z)
	L::Real                   	# half-domain length
	
	# domain definition
	dz::Float64 				# domain spacing
	z::Vector{Float64} 			# domain vector of values (2N + 2 points)

	# magnetic constants 
	B::Float64 					# Bond number
	b::Float64 					# inner rod radius
	E::Float64                  # Bernoulli constant
	
	
	function ferrofluidConstants(N::Int64, L::Real, B::Float64, b::Float64)
		dz = 2*L / (2*N+1)
		z = collect(-Float64(L):dz:Float64(L))
		
		E = 1 - B/2

		new(N, L, dz, z, B, b, E)
	end
end

function c0(k, constants::ferrofluidConstants)
	## linearized wave speed for small amplitude waves c(k)

	B = constants.B
	b = constants.b
	
	c0 = real.(sqrt.(Complex.((1 ./ k).*((-β(1,k,b,1) ./ β(0,k,b,1)) .* (k.^2 .- 1 .+ B)))))

	return c0
	
end

function wall_model(constants::ferrofluidConstants, c::Float64, S::Vector{Float64}, Sz::Vector{Float64}, Szz::Vector{Float64}, Szsq::Vector{Float64})

	## wall model for ferrofluid

	B = constants.B
	E = constants.E

	w = -1 ./ (S.*sqrt.(Szsq)) .+ Szz./(Szsq.^(3/2)) .+ B./(2 .* S.^2) .+ E
	
	return w
	
end

## FU & IL'ICHEV

struct fuConstants <: Constants
	N::Int64  					# number of modes for solution S(z)
	L::Real                   	# half-domain length
	
	# domain definition
	dz::Float64 				# domain spacing
    z::Vector{Float64} 			# domain vector of values (2N + 2 points)

	# physical geometry constants 
	b::Float64 					# inner rod radius
	v::Float64					# fluid velocity 
	
	
	function fuConstants(N::Int64, L::Real, b::Float64, v::Float64)
        dz = 2*L / (2*N+1)
        z = collect(-Float64(L):dz:Float64(L))

        new(N, L, dz, z, b, v)
    end
end

function c0(k, constants::fuConstants)
	## linearized wave speed for small amplitude waves c(k)

	b = constants.b
	v = constants.v

	β0 = besseli(1, k) * besselk(1, k*b) - besseli(1, k*b) * besselk(1, k)
	β1 = besseli(1, k*b) * besselk(0, k) + besseli(0, k) * besselk(1, k*b) - (1/k)*β0

	alpha = (2*β0) / (β0 + k*β1)

    c0 = v * (sqrt(alpha) / (1 + sqrt(alpha)))

	return real.(c0)
	
end

# function wall_model(constants::fuConstants, c::Float64, S::Vector{Float64})

#     v = constants.v
    
#     # Add a small epsilon to prevent division by zero
#     @. w = (1/2) * (1 .- (1 ./ (S.^4))) .* (c .- v).^2

#     return real.(w)
# end

function wall_model(constants::fuConstants, c::Float64, S::Vector{Float64})
    v = constants.v
    w = 0.5 .* (1 .- (1 ./ (S .^ 4))) .* (c - v)^2
    return w
end

## FU & IL'ICHEV (λ2 = 1)

struct fuSimpleConstants <: Constants
	N::Int64  					# number of modes for solution S(z)
	L::Real                   	# half-domain length
	
	# domain definition
	dz::Float64 				# domain spacing
    z::Vector{Float64} 			# domain vector of values (2N + 2 points)

	# physical geometry constants 
	b::Float64 					# inner rod radius
	λ2::Float64					# principal stretches

	vf::Float64					# fluid velocity 
	
	
	function fuSimpleConstants(N::Int64, L::Real, b::Float64, λ2::Float64, vf::Float64)
        dz = 2*L / (2*N+1)
        z = collect(-Float64(L):dz:Float64(L))

        new(N, L, dz, z, b, λ2, vf)
    end
end

function c0(k, constants::fuSimpleConstants)
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

	return real.(c0)
end

function wall_model(constants::fuSimpleConstants, c, S)

	λ2 = constants.λ2
	vf = constants.vf

	w = (λ2^2)/2 .* (1 .- (1 ./ (S.^4))) .* (c - vf/λ2)^2

	return w
end