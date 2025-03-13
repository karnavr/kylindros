using Plots, LaTeXStrings

using SpecialFunctions
using Trapz
using LinearAlgebra
using Statistics

using DelimitedFiles, JSON
using JLD2

using Dates

# include core function files

include("core/models.jl")
include("core/helpers.jl")
include("core/solvers.jl")
include("core/equations.jl")
include("core/bifurcation.jl")
include("core/visualization.jl")
include("core/data.jl")