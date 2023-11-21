#%%
using Pkg
using Revise
Pkg.activate(".")
using FITSIO
using Statistics
#%%
# Read the FITS file
fname = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745/gtexpcube2.fits"

file = FITS(fname, "r")

read(file["ENERGIES"], "Energy")
read(file["ENERGIES"], "Energy")
# get the valeus for energy
earr = Vector{Vector{Float64}}(undef, 10)

for i in 1:10
    earr[i] = read(file["HPXEXPOSURES"], "ENERGY$i")
end

emat = reduce(hcat, earr)

mean(emat)
#%%
fname = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745/gtbin.fits"

file = FITS(fname, "r")

read(file["EBOUNDS"],"E_MIN")
read(file["EBOUNDS"],"E_MAX")

m1 = read(file["SKYMAP"],"CHANNEL1")

sort(m1)