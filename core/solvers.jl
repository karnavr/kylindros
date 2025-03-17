"Custom implementation of numerical solvers for systems of equations"

## WRAPPER FUNCTION

function mySolver(f, initial_guess; solver = :NewtonRaphson, tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	## Solves the system of equations f(x) = 0 using the specified solver (wrapper function)

	if solver == :Newton
		Newton(f, initial_guess, tol = tol, max_iter = max_iter)

	elseif solver == :NewtonRaphson 
	    NewtonRaphson(f, initial_guess, tol = tol, max_iter = max_iter)

	elseif solver == :Broyden
		Broyden(f, initial_guess, tol = tol, max_iter = max_iter)

	elseif solver == :LevenbergMarquardt
		LevenbergMarquardt(f, initial_guess, tol = tol, max_iter = max_iter)

	elseif solver == :Secant
		Secant(f, initial_guess, tol = tol, max_iter = max_iter)

	elseif solver == :JuliaSolver
		

	else
		error("Enter which algorithm you want to use!")
	end
end

## NUMERICAL SOLVERS

function Newton(f, x; tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	for i in 1:max_iter
		J = finite_diff_jacobian(f, x)
		δx = -J \ f(x)  # Newton's update step
		x += δx
		if norm(δx) < tol  # Check for convergence
			return x, i, 0
		end
	end
	
	println("Failed to converge after $max_iter iterations, returning last computed solution.")
	return x, max_iter, 1

end

function NewtonRaphson(f, x; tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	# alpha = 1.0  # Initial step size
	c = 1e-4  # Sufficient decrease parameter
	rho = 0.5  # Step size reduction factor
	
	for i in 1:max_iter
		J = finite_diff_jacobian(f, x)
		δx = -J \ f(x)  # Newton's update step
		t = 1.0  # Initialize step size

		# Backtracking line search
		while norm(f(x + t * δx)) > norm(f(x)) + c * t * dot(f(x), δx)
			t *= rho
		end

		x += t * δx  # Update with the found step size
		if norm(δx) < tol  # Check for convergence
			return x, i, 0
		end
	end
	
	println("Failed to converge after $max_iter iterations, returning last computed solution.")
	return x, max_iter, 1

end

function Broyden(f, x::Vector{Float64}; tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	J = finite_diff_jacobian(f, x)

		for i in 1:max_iter
			δx = -J \ f(x)  # Broyden's update step
			x_new = x + δx
			if norm(δx) < tol  # Check for convergence
				return x, i, 0
			end
			Δf = f(x_new) - f(x)
			J += (Δf - J * δx) * δx' / (δx' * δx)  # Update Jacobian approximation
			x = x_new
		end

		println("Failed to converge after $max_iter iterations, returning last computed solution.")
		return x, max_iter, 1

end

function LevenbergMarquardt(f, x::Vector{Float64}; tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	λ = 1e-3  # Initial damping factor
		ν = 2.0   # Damping factor increment

		for i in 1:max_iter
			J = finite_diff_jacobian(f, x)
			A = J' * J + λ * I  # Levenberg-Marquardt update
			g = J' * f(x)
			δx = -A \ g

			x_new = x + δx
			if norm(f(x_new)) < norm(f(x))  # If the step improves the solution
				x = x_new
				λ /= ν
			else
				λ *= ν
			end

			if norm(δx) < tol  # Check for convergence
				return x, i, 0
			end
		end
		
		println("Failed to converge after $max_iter iterations, returning last computed solution.")
		return x, max_iter, 1

end

function Secant(f, x::Vector{Float64}; tol::Float64 = 1e-8, max_iter::Int64 = 1000)

	δx = ones(length(x)) * 1e-4
		for i in 1:max_iter
			f_val = f(x)
			if norm(f_val) < tol
				return x
			end

			J_approx = finite_diff_jacobian(f, x)
			δx = -J_approx \ f_val
			x += δx

			if norm(δx) < tol  # Check for convergence
				return x, i
			end
		end
		error("Failed to converge after $max_iter iterations")

end

## AUXILIARY FUNCTIONS

function finite_diff_jacobian(f, x)

	## Computes the Jacobian of the function f at the point x using finite differences

    h = 1e-8  # Small perturbation (typically sqrt(machine epsilon))
    n = length(x)
    J = zeros(n, n)
    fx = f(x)

    Threads.@threads for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += h
        J[:, i] = (f(x_perturbed) - fx) / h
    end
    return J
end

function better_finite_diff_jacobian(f, x)

	# computes the jacobian of f at x using finite differences from the 
	# actual finite difference package 

	return FiniteDiff.finite_difference_jacobian(f, x)

end

function auto_diff_jacobian(f, x)

	## Computes the Jacobian of the function f at the point x using automatic differentiation

	return ForwardDiff.jacobian(f, x)

end

# Helper function to safely convert Complex tracked values
function safe_complex(z)
    if z isa Complex
        return Complex(real(z), imag(z))
    else
        return z
    end
end

# Simpler approach for Bessel functions with tracked values
function SpecialFunctions.besseli(nu::Real, z::Complex{<:ReverseDiff.TrackedReal})
    # Extract real and imaginary parts of z
    re_z = real(z)
    im_z = imag(z)
    
    # Convert to regular values
    z_val = Complex(ReverseDiff.value(re_z), ReverseDiff.value(im_z))
    
    # Compute function value
    result = besseli(nu, z_val)
    
    # Just return the result without trying to track the complex directly
    # This will allow ReverseDiff to handle the derivative propagation internally
    return result
end

function SpecialFunctions.besselk(nu::Real, z::Complex{<:ReverseDiff.TrackedReal})
    # Extract real and imaginary parts of z
    re_z = real(z)
    im_z = imag(z)
    
    # Convert to regular values
    z_val = Complex(ReverseDiff.value(re_z), ReverseDiff.value(im_z))
    
    # Compute function value
    result = besselk(nu, z_val)
    
    # Just return the result without trying to track the complex directly
    return result
end