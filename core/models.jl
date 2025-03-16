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
	
	c0 = sqrt.((1 ./ k).*((-β(1,k,b,1) ./ β(0,k,b,1)) .* (k.^2 .- 1 .+ B)))

	return c0
	
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
	λ1::Float64					# principal stretch 1
	λ2::Float64					# principal stretch 2
	
	vf::Float64					# fluid velocity 
	
	
	function fuConstants(N::Int64, L::Real, b::Float64, λ1::Float64, λ2::Float64, vf::Float64)
        dz = 2*L / (2*N+1)
        z = collect(-Float64(L):dz:Float64(L))

        new(N, L, dz, z, b, λ1, λ2, vf)
    end
end

function c0(k, constants::fuConstants)
	## linearized wave speed for small amplitude waves c(k)

	b = constants.b
	λ1 = constants.λ1
	λ2 = constants.λ2
	vf = constants.vf

	β0 = besseli(1, k) * besselk(1, k*b) - besseli(1, k*b) * besselk(1, k)
	β1 = besseli(1, k*b) * besselk(0, k) + besseli(0, k) * besselk(1, k*b) - (1/k)*β0

	# quadratic coeffs (Ac^2 + Bc + D) 
	A = k*β1 - ((2*λ2^2 * β0) ./ λ1) + β0
	B = (4 * vf * λ2 * β0) ./ λ1
	D = -2 * vf^2 * β0 ./ λ1

	# solve quadratic equation 
	c0 = solve_quadratic(A, B, D)

	return real(c0)
	
end

function wall_model(constants::fuConstants, c::Float64, S::Vector{Float64})

	λ2 = constants.λ2
	vf = constants.vf

	w = (λ2^2)/2 .* (1 .- (((λ1^4) ./ ((λ1 .+ S .- 1).^4))) .* (c - vf/λ2)^2)

	return w
end

function unpackConstants(constants::fuConstants)

	# unpacks values and returns as string 
	N = constants.N
	L = constants.L
	b = constants.b
	λ1 = constants.λ1
	λ2 = constants.λ2
	vf = constants.vf
	
	# create a string of all the values 
	return string(N, L, b, λ1, λ2, vf)

end

## FU & IL'ICHEV (λ1 = 1)

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

	return real(c0)
	
end

function wall_model(constants::fuSimpleConstants, c, S::Vector)

	λ2 = constants.λ2
	vf = constants.vf

	w = (λ2^2)/2 .* (1 .- (1 ./ (S.^4))) .* (c - vf/λ2)^2

	return w
end

function unpackConstants(constants::fuSimpleConstants)

	# unpacks values and returns as string 
	N = constants.N
	L = constants.L
	b = constants.b
	λ2 = constants.λ2
	vf = constants.vf
	
	# create a string of all the values 
	return string(N, L, b, λ2, vf)

end