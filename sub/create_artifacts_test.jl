using Pkg
using Revise
Pkg.activate(".")

using MapGenData
using HDF5
using Unitful
using StaticArrays
using JLD2
using Base.Threads
using Healpix
using ProgressMeter
using QuadGK
using BenchmarkTools
using MultiQuad
using Logging

using PyCall
hp = pyimport("healpy")

using PyPlot
#%%
# MapGenData.clear_cache()
#%%
@info "Using nthreads = $(nthreads())"

hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_0.5_1000_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts"
artifact_cache = MapGenData.artifact_cache

fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

@info "Start creating FITS artifacts"
make_fits_artifact(fits_artifact)

Earr = [500, 1000, 2000, 5000, 10_000, 50_000, 200_000, 1_000_000] * u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])

#%%
jld2_artifact = JLD2Artifact("./", 1024, Emin_macro, Emax_macro)
fg5 = MapGenData.read_galactic_fg_v05(jld2_artifact)

#%%
plt.clf()
hp.mollview(fg5[1][:,end], title="Galactic foreground", unit="flux", norm="log")
plt.gcf()
# MapGenData.write_gf_v07_map_smoothed_as_jld2(jld2_artifact)
#%%
# MapGenData.write_gf_v07_counts_map_as_jld2("./", jld2_artifact; compress=true)