include("../functions.jl")
println("Avaliable threads: ", Threads.nthreads())

T = ArbReal
setprecision(T, 52); setextrabits(10);

## Initialize 

branchN = 200
a1Vals = collect(range(T(0.01), T(0.29), branchN + 1))

# Define the values for the parameters
N = 36
L = T(π)
b = T(0.1)
λ2 = T(1.5)
vf = T(1.0)

# Create an instance of the Constants struct
constants = fuSimpleConstants(T, N, L, b, λ2, vf)

# initialize wave speed and wave number 
k1 = T(1)*T(π)/constants.L
cInitial = c0(k1, constants)

# initial guess 
initial_guess = T(1e-16).*ones(T, branchN+1, constants.N+2)
initial_guess[1,1:4] = T.([cInitial, 1.0, a1Vals[1], 1e-16])

## Compute solution branch
@time solutions, constants, metadata = bifurcation(initial_guess, a1Vals, branchN, constants, tol = T(1e-12), solver = :NLSolver, max_iter = 10000, overwrite = true);