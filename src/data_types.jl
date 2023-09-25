struct FITSArtifact
    fitsdir::String
    outdir::String
end

struct JLD2Artifact{S}
    outdir::String
    nside::Int
    Emin_array::SVector{S, Float64} # MeV
    Emax_array::SVector{S, Float64} # MeV
end

function JLD2Artifact{S}(outdir::String, nside::Int, Emin_array::AbstractVector{Float64}, Emax_array::AbstractVector{Float64}) where {S}
    @assert nside <= 2048 "nside must be <= 2048"
    @assert ispow2(nside) "nside must be a power of 2"
    Emin_sv = SVector{S, Float64}(Emin_array)
    Emax_sv = SVector{S, Float64}(Emax_array)

    # check wheter macrobins are consistent with microbins
    Emin_micro, Emax_micro = get_E_bins()
    @assert check_bins(Emin_micro, Emax_micro, Emin_sv, Emax_sv) "Energy bins are not consistent, make sure that the macrobins can be interpreted as a combination of the microbins contained in the fits files"
    
    return JLD2Artifact{S}(outdir, nside, Emin_sv, Emax_sv)
end

function JLD2Artifact(outdir::String, nside::Int, Emin_array::AbstractVector, Emax_array::AbstractVector)
    return JLD2Artifact{length(Emin_array)}(outdir, nside, Emin_array, Emax_array)
end

function JLD2Artifact(outdir::String, nside::Int, Emin_array::AbstractVector{T}, Emax_array::AbstractVector{T}) where {T<:Energy}
    Emin_ = ustrip.(u"MeV", Emin_array)
    Emax_ = ustrip.(u"MeV", Emax_array)
    return JLD2Artifact{length(Emin_array)}(outdir, nside, Emin_, Emax_)
end
