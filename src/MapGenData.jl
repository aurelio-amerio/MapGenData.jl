module MapGenData

export make_fits_artifact, make_jld2_artifacts, FITSArtifact, JLD2Artifact

artifact_cache = ""

include("includes.jl")
include("utils.jl")
include("data_types.jl")
include("process_fits.jl")
include("process_galactic_foreground.jl")
include("make_artifacts.jl")



function __init__()
    global artifact_cache = @get_scratch!("artifact_cache")
end

end
