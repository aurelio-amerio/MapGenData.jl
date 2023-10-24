using Pkg
using Revise
Pkg.activate(".")

ENV["MapGenData_cache_label"] = "1bin"

using MapGenData
using HDF5
using Unitful
using StaticArrays
using JLD2
using Base.Threads
#%%
MapGenData.clear_cache()
#%%
@info "Using nthreads = $(nthreads())"

hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_0.5_1000_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts_1bin"
artifact_cache = MapGenData.artifact_cache

fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

@info "Start creating FITS artifacts"
make_fits_artifact(fits_artifact)

Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])

#%%
@info "Start creating jld artifacts"
# first compute the foreground at nside 1024
let
    jld2_artifact = JLD2Artifact("./", 1024, Emin_macro, Emax_macro)

    MapGenData.write_gf_map_smoothed_as_jld2(jld2_artifact, version=7)
    MapGenData.write_gf_map_smoothed_as_jld2(jld2_artifact, version=5)
end


for nside in [1024, 64, 128, 256, 512]
    jld2_artifact = JLD2Artifact(artifacts_folder, nside, Emin_macro, Emax_macro)
    make_jld2_artifacts(jld2_artifact_)
end


