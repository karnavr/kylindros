

function createResultsDict(constants::Constants, solutions::Matrix{Float64}, metadata::Dict{String, Any})

    results_dict = Dict(
        "solutions" => solutions,
        "constants" => constants,
        "metadata" => metadata
    )

    return results_dict
end

function generateFileName(results_dict::Dict{String, Any})

    ## creates a unique filename based on the model name, solver type, tolerance, and everything in constants 

    model_name = results_dict["metadata"]["model"]
    solver_type = results_dict["metadata"]["solver"]

    constants = results_dict["constants"]

    N = constants.N
    tol = results_dict["metadata"]["tol"]
    branchN = results_dict["metadata"]["branchN"]
    

    # create a unique filename based on the model name, solver type, tolerance, and everything in constants 
    filename = "$(model_name)_$(N)_$(tol)_$(branchN)_$(solver_type)"

    return filename
end 

function getModelName(constants::Constants)
    
    # Get the full type name as string and remove "Constants" from the end
    full_name = string(nameof(typeof(constants)))

    model_name = full_name[1:end-9] # Remove "Constants" (9 characters)

    return model_name
end