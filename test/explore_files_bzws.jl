#%%
using Pkg
using Revise
Pkg.activate(".")
using FITSIO
using Statistics
#%%
# Read the FITS file gtmktime
fname = "/data/users/Aurelio/fermi/output/sourceveto_nside2048_front_1_10GeV/gtmktime.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")
gti_start[79473]

file[1]
#%% gtselect
fname = "/data/users/Aurelio/fermi/output/sourceveto_nside2048_front_1_10GeV/gtselect.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")
gti_start[end]

file[1]
# %% gtbin
fname = "/data/users/Aurelio/fermi/output/sourceveto_nside2048_front_1_10GeV/gtbin.fits"
file = FITS(fname, "r")

file

file["GTI"]

gti_start=read(file["GTI"], "START")

c1=read(file["SKYMAP"],"CHANNEL1")