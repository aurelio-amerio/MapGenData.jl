#%%
using Pkg
using Revise
Pkg.activate(".")
using FITSIO
using Statistics
#%% gtselect
fname = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745_v2/gtselect.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")

gti_start[end]

file[1]
#%%
# Read the FITS file gtmktime
fname = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745_v2/gtmktime.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")

gti_start[79476]

file[1]
#%% gtbin
fname = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745_v2/gtbin.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")

gti_start[79476]

file[1]

c1=read(file["SKYMAP"],"CHANNEL1")