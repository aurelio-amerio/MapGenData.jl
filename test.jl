using Pkg
using Revise
Pkg.activate(".")

using MapGenData

#%%
hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_200_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts"

fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

# mkpath(fits_artifact.outdir)

# MapGenData.make_fits_artifact(fits_artifact)
#%%
# MapGenData.delete_cache()


#%%

#%%


using HDF5


Emin_micro, Emax_micro = MapGenData.get_E_bins()


Emin_micro
Emax_micro

Emin_macro = Emin_micro[1:5:end]
Emax_macro = Emax_micro[5:5:end]
jld2_artifact = JLD2Artifact(artifacts_folder, 128, "test1", Emin_macro, Emax_macro)

#%%
map_fermi = MapGenData.read_fermi_map(jld2_artifact)
pwd()
MapGenData.write_fermi_map_as_jld2("$(pwd())", jld2_artifact)
#%%
using Unitful
MapGenData.get_E_bins()
MapGenData.get_E_bins(u"GeV")
#%%
filepath = joinpath(MapGenData.artifact_cache, "fits", "gtexpcube2.h5")

fid = h5open(filepath, "r")

keys(fid["HPXEXPOSURES/DATA"])
close(fid)


MapGenData.write_exposure_map_as_jld2("$(pwd())", jld2_artifact)