using Pkg
using Revise
Pkg.activate(".")

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
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts"
artifact_cache = MapGenData.artifact_cache

fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

make_fits_artifact(fits_artifact)

Earr = [500, 1000, 2000, 5000, 10_000, 50_000, 200_000, 1_000_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])

#%%
for nside in [1024, 64, 128, 256, 512]
    jld2_artifact = JLD2Artifact(artifacts_folder, nside, Emin_macro, Emax_macro)
    make_jld2_artifacts(jld2_artifact)
end