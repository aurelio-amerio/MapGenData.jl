struct FITSArtifact
    fitsdir::String
    outdir::String
end

struct JLD2Artifact
    fitsdir::String
    outdir::String
    nside::Int
    name::String
    Emin::typeof(1.0u"GeV")
    Emax::typeof(1.0u"GeV")
    nmicrobins::Int
    nbins::Int
end