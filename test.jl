using Pkg
using Revise
Pkg.activate(".")

using MapGenData

#%%
fits_artifact = FITSArtifact("/home/aure/Desktop/data", "/home/aure/Desktop/data/out")

mkpath(fits_artifact.outdir)

MapGenData.make_fits_artifact(fits_artifact)
#%%
using ArtifactUtils
tmp_dir = "/tmp/tmp_art"
artifact_from_directory(tmp_dir)

using Downloads
gll_iem_url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits"
Downloads.download(gll_iem_url, joinpath(tmp_dir, "gll_iem_v07.fits"))