"Functions for solution data saving and loading"

function generateFileName(metadata::Dict{String, Any})

    "creates a unique filename based on the model name and current date/time"

    # get the model name (everything before "Constants" in the full name)
    model_name = metadata["model"]

    date_string = metadata["timestamp"]

    # create a unique filename based on the model name  
    filename = "$(model_name)_$(date_string)"

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

    # format the date and time as a string (day-month_hour-minute-second_year)
    date_string = Dates.format(current_time, "dd-mm_HH-MM-SS_yyyy")

    return date_string
end

function solutionExists(constants::Constants, branchN::Int64, tol::Float64) 

    "Loops through all the files in the results directory and checks if a solution exists for the given constants. 
    
    If it does, it returns the filename. If it doesn't, it returns false."

    # get the model name
    model_name = getModelName(constants)
    
    # Check if the directory exists
    if !isdir("results/$(model_name)")
        return false
    end

    # loop through all the files in the results/model_name directory
    for file in readdir("results/$(model_name)")

        if endswith(file, ".jld2")
            
            # load the file
            saved_results = load("results/$(model_name)/$(file)")
            
            # Now we can directly access the components
            if haskey(saved_results, "constants") && haskey(saved_results, "metadata")

                # check if the constants are the same
                if compareConstants(constants, saved_results["constants"])

                    # check if the tolerance and branchN are the same
                    if saved_results["metadata"]["tol"] == tol && 
                       saved_results["metadata"]["branchN"] == branchN
                        return file
                    end

                end

            end
            
        end
    end

    return false
end

function compareConstants(constants1::Constants, constants2::Constants)

    "Compares two constants structs for all fields except for dz and the z vector and returns true if they are equal, false otherwise."

    # get the fieldnames of the Constants struct
    all_fieldnames = fieldnames(typeof(constants1))

    # loop through all the fieldnames
    for field in all_fieldnames
        if field != :dz && field != :z
            if getfield(constants1, field) != getfield(constants2, field)
                return false
            end
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