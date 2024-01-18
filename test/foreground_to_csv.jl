using Pkg
using Revise
Pkg.activate(".")
using MapGenData
using FITSIO
using Unitful
using Healpix
using DataFrames
using CSV
using ProgressMeter
#%%
filepath = joinpath(MapGenData.fits_cache, "gll_iem_v07.fits")
file = FITS(filepath, "r")
energy_fg1 = read(file[2], "energy")
model_ = Float64.(read(file[1]))
close(file)

model = permutedims(model_, (3,2,1))
eneb_fg = length(energy_fg1)

#%%
nside=512
res = Resolution(nside)
Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])

jld2_artifact = JLD2Artifact("./", nside, Emin_macro, Emax_macro)


model, en = MapGenData.read_galactic_fg_v07(jld2_artifact)
model
npix = 12*nside^2
pixels = 1:npix
theta = zeros(npix)
phi = zeros(npix)
for i in pixels
    theta[i], phi[i] = pix2angRing(res, i)
end

en
lat2colat.(deg2rad.(glat)), deg2rad.(glon) # theta, phi

glat = rad2deg.(colat2lat.(theta))
glon = rad2deg.(phi)

# cm^-2 s^-1 sr^-1 MeV^-1
# create a dataframe with the following column names: 
# gLat, gLon, and 28 columns, one for each energy bin
#%%
df = DataFrame("gLat [deg]"=>glat, "gLon [deg]"=>glon)

CSV.write("gll_iem_v7_csv/lonlat.csv", df)

df = DataFrame("energies [MeV]" => en)

CSV.write("gll_iem_v7_csv/energies.csv", df)

@showprogress for i in eachindex(en)
    df = DataFrame("foreground [cm^-2 s^-1 sr^-1 MeV^-1]" => model[:,i])

    CSV.write("gll_iem_v7_csv/gll_iem_v7_$(lpad(i, 2, '0')).csv", df)
end
#%%
