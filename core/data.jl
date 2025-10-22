"Functions for solution data saving and loading"

using Logging

function results_dir(subdir::String = "")

    "Return the absolute path to the results directory, optionally within a subdirectory."

    base = joinpath(dirname(@__DIR__), "results")
    if isempty(subdir)
        return base
    else
        return joinpath(base, subdir)
    end
end

function generateFileName(metadata::Dict{String, Any}, constants::Constants)

    "creates a unique filename based on the model name, date/time, and number of Fourier modes"

    # extract components
    model_name = metadata["model"]
    date_string = metadata["timestamp"]

    # create filename: model_timestamp_N##
    filename = "$(model_name)_$(date_string)_N$(constants.N)"

    return filename
end 

function getModelName(constants::Constants)
    
    # Get the full type name as string and remove "Constants" from the end
    full_name = string(nameof(typeof(constants)))

    model_name = full_name[1:end-9] # Remove "Constants" (9 characters)

    return model_name
end

function getCurrentDateToString()
    # get the current date and time
    current_time = Dates.now()

    # format the date and time as a compact string (ddmmHHMM)
    date_string = Dates.format(current_time, "ddmmHHMM")

    return date_string
end

function getSaveDirectory(constants::Constants, save_dir::String)
    
    "Determines the save directory: auto mode (adds model subdir under library/) or custom mode (use as-is)."
    
    model_name = getModelName(constants)
    return isempty(save_dir) ? joinpath(results_dir("library"), model_name) : results_dir(save_dir), model_name
end


function compareConstants(constants1, constants2)

    "Compares two constants structs for all fields except for dz and the z vector and returns true if they are equal, false otherwise."

    # get the fieldnames of both structs
    fieldnames1 = fieldnames(typeof(constants1))
    fieldnames2 = fieldnames(typeof(constants2))
    
    # check if they have the same fields (excluding dz and z)
    excluded = [:dz, :z]
    fields1 = filter(f -> f ∉ excluded, fieldnames1)
    fields2 = filter(f -> f ∉ excluded, fieldnames2)
    
    # if different fields, not equal
    fields1 != fields2 && return false

    # loop through all the comparable fieldnames
    for field in fields1
        if getfield(constants1, field) != getfield(constants2, field)
            return false
        end
    end

    return true
end

function readSolution(file_path::String)

    "Reads a solution from a file and returns the constants and metadata."

    # load the file
    saved_results = load(file_path)
    
    return saved_results["solutions"], saved_results["constants"], saved_results["metadata"]
end

function findMatchingSolutions(dir::String, constants::Constants, branchN::Int64, a1Vals::Vector{Float64}, tol::Float64, solver::Symbol)

    "Scans directory and returns all filenames matching the given parameters."

    matching_files = String[]
    
    for file in readdir(dir)
        !endswith(file, ".jld2") && continue
        
        saved_results = Logging.with_logger(Logging.NullLogger()) do
            load(joinpath(dir, file))
        end
        
        (!haskey(saved_results, "constants") || !haskey(saved_results, "metadata")) && continue
        compareConstants(constants, saved_results["constants"]) || continue
        
        meta = saved_results["metadata"]
        if meta["tol"] == tol && 
           meta["branchN"] == branchN &&
           all(isapprox.(extrema(a1Vals), extrema(meta["a1Vals"]))) &&
           meta["solver"] == solver
            push!(matching_files, file)
        end
    end

    return matching_files
end

function handleExistingSolution(constants::Constants, branchN::Int64, a1Vals::Vector{Float64}, tol::Float64, solver::Symbol, save_dir::String, overwrite::Bool)

    "Handles existing solution loading or deletion based on overwrite flag. Returns (existing_solution, save_directory, model_name)."

    # determine save directory and model name
    dir, model_name = getSaveDirectory(constants, save_dir)

    # check if directory exists
    if !isdir(dir)
        return (nothing, dir, model_name)
    end

    # find all matching files
    matching_files = findMatchingSolutions(dir, constants, branchN, a1Vals, tol, solver)

    # handle matches based on overwrite flag
    if !isempty(matching_files)
        if !overwrite
            # load and return first match
            full_path = joinpath(dir, matching_files[1])
            println("Solution branch already exists.")
            println("File: $(full_path)")
            return (readSolution(full_path), dir, model_name)
        else
            # delete all matches
            println("Found $(length(matching_files)) duplicate saved results. Deleting (overwrite = true) and re-computing.")
            for file in matching_files
                rm(joinpath(dir, file))
            end
        end
    end

    return (nothing, dir, model_name)
end