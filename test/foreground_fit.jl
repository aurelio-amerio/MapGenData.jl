using Pkg
using Revise
Pkg.activate(".")

ENV["MapGenData_cache_label"] = "1bin_1_10GeV"

using MapGenData

using Base.Threads
using HDF5
using Healpix
using FITSIO
using Interpolations
using JLD2
using ProgressMeter
using PyCall
using QuadGK
using StaticArrays
using Unitful
using Unitful: Energy

Minuit = pyimport("iminuit").Minuit
#%%
function pl_f(Agal::Real, Fiso::Real, fore_map::AbstractArray, k_arr::AbstractArray, fac::AbstractArray)
    lambda_arr = Agal*fore_map.* fac .+ Fiso* fac # valore previsto 
    llarg = -k_arr .* log.(lambda_arr) .+ lambda_arr
    return sum(llarg[.!isnan.(llarg)])
end

function fit_foreground(;Agal_0=1.0, Fiso_0=1e-7, limit_Agal=[0.2,2.0], limit_Fiso=[1e-15,1e-6], fermi_map, GF, exposure_map, srpixel, pls_mask) where {T}
    @assert length(fermi_map) == length(GF) == length(exposure_map) "shape mismatch between the provided maps"
    
    fac = view(Matrix(exposure_map .* srpixel), pls_mask,:)
    fore_map = view(Matrix(GF), pls_mask, :)
    k_arr = view(Matrix(fermi_map), pls_mask, :)


    nchannels = size(k_arr)[2]

    Agal_vals = zeros(nchannels)
    Fiso_vals = zeros(nchannels)
    Agal_errors = zeros(nchannels)
    Fiso_errors = zeros(nchannels)

    function _f(Agal, Fiso, bin) 
        return pl_f(Agal, Fiso, fore_map[:,bin], k_arr[:,bin], fac[:,bin])
    end

    for i in 1:nchannels
        m = Minuit((Agal,Fiso)->_f(Agal,Fiso,i), Agal_0, Fiso_0)
        m.limits[1] = limit_Agal
        m.limits[2] = limit_Fiso
        m.errordef = Minuit.LIKELIHOOD
        m.migrad()
        m.hesse()

        vals = Float64.(collect(m.values))
        errors = Float64.(collect(m.errors))

        Agal_vals[i] = vals[1]
        Fiso_vals[i] = vals[2]
        Agal_errors[i] = errors[1]
        Fiso_errors[i] = errors[2]
    end

    return (Agal_vals, Agal_errors), (Fiso_vals, Fiso_errors)
end

function read_galactic_fg_v07(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    res = Resolution(nside)
    npix = 12 * nside^2

    # galactic diffuse map
    filepath = joinpath(MapGenData.fits_cache, "gll_iem_v07.fits")

    file = FITS(filepath, "r")
    energy_fg1 = read(file[2], "energy")
    model_ = Float64.(read(file[1]))
    close(file)

    model = permutedims(model_, (3,2,1))
    eneb_fg = length(energy_fg1)

    dec2 = zeros(npix)
    ra2 = zeros(npix)
    for j in 1:npix
        (dec2[j], ra2[j]) = pix2angRing(res, j)
    end
    nres = 8 # 8
    dec2 = ( -dec2 * 180 / pi .+ 180 ) .* nres
    ra2 = ( ( -ra2 * 180 / pi .+ 360 .+ 180 ) .% 360 ) .* nres

    model_heal = ones((npix,eneb_fg))

    ra2[ra2.>2879].=2879
    dec2[dec2.>1440].=1440

    for ii in eachindex(energy_fg1)
        model_heal[:,ii] .= MapGenData.map_coordinates(model[ii,:,:],[dec2,ra2]) 
    end

    return model_heal, energy_fg1
end

function process_galactic_fg_v07(jld2_artifact, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    
    model_heal, energy_fg1 = read_galactic_fg_v07(jld2_artifact)

    Emin = ustrip(u"MeV", Emin) 
    Emax = ustrip(u"MeV", Emax) 
    npix = 12*nside^2
    map_gf = HealpixMap{Float64, RingOrder}(nside)
    p = Progress(npix)
    @threads for j in 1:npix
        interfunc_tmp = linear_interpolation(log10.(energy_fg1), log10.(model_heal[j,:]))

        interfunc(E) = 10 ^ interfunc_tmp(log10(E))

        map_gf[j] = quadgk(interfunc, Emin, Emax, rtol=1e-3)[1]
        next!(p)
    end
    return map_gf
end

function read_galactic_fg_v05(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    res = Resolution(nside)
    npix = 12 * nside^2

    # galactic diffuse map
    filepath = joinpath(MapGenData.fits_cache, "gll_iem_v05_rev1.fit")
    file = FITS(filepath, "r")
    energy_fg1 = read(file[2], "Energy")
    model_ = Float64.(read(file[1]))
    close(file)

    model = permutedims(model_, (3,2,1))
    eneb_fg = length(energy_fg1)

    dec2 = zeros(npix)
    ra2 = zeros(npix)
    for j in 1:npix
        (dec2[j], ra2[j]) = pix2angRing(res, j)
    end
    # nres = Int(log2(nside))#8 # 8
    nres = 8 # 8
    dec2 = ( -dec2 * 180 / pi .+ 180 ) .* nres
    ra2 = ( ( +ra2 * 180 / pi .+ 360 .+ 180 ) .% 360 ) .* nres

    model_heal = ones((npix,eneb_fg))

    ra2[ra2.>2879].=2879
    dec2[dec2.>1440].=1440

    for ii in eachindex(energy_fg1)
        model_heal[:,ii] .= MapGenData.map_coordinates(model[ii,:,:],[dec2,ra2]) 
    end

    return model_heal, energy_fg1
end

function process_galactic_fg_v05(jld2_artifact, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    
    model_heal, energy_fg1 = read_galactic_fg_v05(jld2_artifact)

    Emin = ustrip(u"MeV", Emin) 
    Emax = ustrip(u"MeV", Emax) 
    npix = 12*nside^2
    map_gf = HealpixMap{Float64, RingOrder}(nside)
    p = Progress(npix)
    @threads for j in 1:npix
        interfunc_tmp = linear_interpolation(log10.(energy_fg1), log10.(model_heal[j,:]))

        interfunc(E) = 10 ^ interfunc_tmp(log10(E))

        map_gf[j] = quadgk(interfunc, Emin, Emax, rtol=1e-3)[1]
        next!(p)
    end
    return map_gf
end
#%%
# MapGenData.clear_cache()
#%%
@info "Using nthreads = $(nthreads())"

# hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_0.5_1000_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts_test"
artifact_cache = MapGenData.artifact_cache

# fits_artifact = FITSArtifact(hdf5_folder, artifacts_folder)

# @info "Start creating FITS artifacts"
# MapGenData.fetch_fermilat_data(fits_artifact)
# MapGenData.make_fits_artifact(fits_artifact)

Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])
#%%
fermi_map = load(joinpath(artifact_cache,"nside1024","fermi_map.jld2"))["fermi_map"]
exposure_map = load(joinpath(artifact_cache,"nside1024","exposure_map.jld2"))["exposure_map"]
pls_mask = load("test/pls_mask_1024_2deg.jld2")["pls_mask"]
#%%
jld2_artifact = JLD2Artifact(artifacts_folder, 1024 , Emin_macro, Emax_macro)
#%%
fg_map_v7 = [process_galactic_fg_v07(jld2_artifact, 1024)]
fg_map_v5 = [process_galactic_fg_v05(jld2_artifact, 1024)]
#%%
fit_foreground(fermi_map=fermi_map, GF=fg_map_v7, exposure_map=exposure_map, srpixel=4pi/(12*1024^2), pls_mask=pls_mask)
fit_foreground(fermi_map=fermi_map, GF=fg_map_v5, exposure_map=exposure_map, srpixel=4pi/(12*1024^2), pls_mask=pls_mask)
#%%
filepath = joinpath(MapGenData.fits_cache, "gll_iem_v05_rev1.fit")
file = FITS(filepath, "r")
energy_fg1 = read(file[2], "Energy")
model_ = Float64.(read(file[1]))
close(file)

model = permutedims(model_, (3,2,1))
sqrt(1441*2880 / 12)
#%%
fermi_map_64 = MapGenData.ud_grade(fermi_map[1], 64, power=-2)
exposure_map_64 = MapGenData.ud_grade(exposure_map[1], 64)
fg_map_v7_64 = MapGenData.ud_grade(fg_map_v7[1], 64)
fg_map_v5_64 = MapGenData.ud_grade(fg_map_v5[1], 64)
pls_mask_64 = load("test/pls_mask_64.jld2")["pls_mask"]
#%%

fit_foreground(fermi_map=[fermi_map_64], GF=[fg_map_v5_64], exposure_map=[exposure_map_64], srpixel=4pi/(12*64^2), pls_mask=pls_mask_64)