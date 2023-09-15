using Pkg
using Revise
Pkg.activate(".")

using MapGenData
using HDF5
using Unitful
using StaticArrays

#%%
hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_200_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts"
artifact_cache = MapGenData.artifact_cache

fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

Emin_micro, Emax_micro = MapGenData.get_E_bins(u"MeV")

Emin_macro = Emin_micro[1:5:end]
Emax_macro = Emax_micro[5:5:end]
jld2_artifact = JLD2Artifact(artifacts_folder, 1024, "test1", Emin_macro, Emax_macro)

#%%
MapGenData.write_gf_v07_map_smoothed_as_jld2(jld2_artifact)

#%%
