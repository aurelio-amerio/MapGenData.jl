#%%
using Pkg
using Revise
Pkg.activate(".")
using MapGenData
using FITSIO
using Statistics
using Unitful
using Healpix
using ProgressMeter
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
#%%
using JLD2
ENV["MapGenData_cache_label"] = "1bin_1_10GeV_xgam"
using MapGenData

MapGenData.artifact_cache
data = load(joinpath(MapGenData.artifact_cache, "nside1024","galactic_foreground_v07_smoothed_counts.jld2" ))

data["gf"]
#%% test galactic foreground
nside=1024
Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])
jld2_artifact = JLD2Artifact("./", nside, Emin_macro, Emax_macro)
model_heal, energy_fg1 = MapGenData.read_galactic_fg_v07(jld2_artifact)

model_heal

i=1
Emin = jld2_artifact.Emin_array[i]*u"MeV"
Emax = jld2_artifact.Emax_array[i]*u"MeV"
lmax=nside*4
PSF_theta = MapGenData.get_PSF_theta(jld2_artifact)
            
PSF_theta_bini(theta) = PSF_theta(theta, i)
model_smoothed, energy_fg1_filtered = MapGenData.convolve_fg_model_with_PSF(model_heal, lmax, PSF_theta_bini, energy_fg1; Emin=Emin, Emax=Emax)

model_smoothed
#%% counts map
# now we interpolate the model in energy and store the interpolations in an array
gf_model_interpolated = Vector{Function}(undef, length(jld2_artifact.Emin_array))


gf_model_interpolated[1], _ = MapGenData.interpolate_galactic_fg(model_smoothed, energy_fg1_filtered, jld2_artifact)


# now we compute the integral of the model in each energy bin in order to obtain the foreground model in units of counts
gf_integral = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))

exp_map_E, _ = MapGenData.get_exposure_map_interpolation(jld2_artifact)
@info "Convolving the smoothed foreground template with the exposure map and computing the integral in each energy bin"
p = Progress(length(gf_integral))
for i in eachindex(gf_integral)
    Emin = jld2_artifact.Emin_array[i]*u"MeV"
    Emax = jld2_artifact.Emax_array[i]*u"MeV"
    gf_integral[i] = MapGenData.process_galactic_fg_smoothed_counts(gf_model_interpolated[i], exp_map_E, jld2_artifact; Emin=Emin, Emax=Emax)
    next!(p)
end


gf_integral[1][1]

