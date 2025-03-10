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