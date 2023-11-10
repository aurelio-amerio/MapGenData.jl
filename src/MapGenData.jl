module MapGenData

export make_fits_artifact, make_jld2_artifacts, FITSArtifact, JLD2Artifact

artifact_cache = ""
artifact_cache_name = ""

fits_cache = ""


include("includes.jl")
include("utils.jl")
include("data_types.jl")
include("process_fits.jl")
include("process_galactic_foreground.jl")
include("make_artifacts.jl")



function __init__()
    if "MapGenData_cache_label" in keys(ENV)
        global artifact_cache_name = "artifact_cache_$(ENV["MapGenData_cache_label"])"
    else
        global artifact_cache_name = "artifact_cache"
    end
    global artifact_cache = @get_scratch!(artifact_cache_name)

    global fits_cache = @get_scratch!("fits_cache")
    return
end

end
