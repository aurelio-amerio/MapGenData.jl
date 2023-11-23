using Pkg
using Revise
Pkg.activate(".")
using MapGenData
using Base.Threads
using BenchmarkTools
using FITSIO
using Healpix
using Interpolations
using LegendrePolynomials
using ProgressMeter 
using QuadGK
using Statistics
using Unitful
using Unitful: Energy
using Memoize
using JLD2

using PyCall
Minuit = pyimport("iminuit").Minuit
#%% using library
function get_exposure_map_interpolation_modified(jld2_artifact::JLD2Artifact)

    nside = jld2_artifact.nside
    filepath = joinpath(MapGenData.fits_cache, "gtexpcube2.h5")

    Emin_micro, Emax_micro = get_E_bins()
    En_arr = (Emin_micro .+ Emax_micro) ./ 2

    map_binned = _read_exposure_map_helper(nside, filepath, Emin_micro, Emax_micro, Emin_micro, Emax_micro) 
    map_for_itp = convert(Matrix{Float64}, map_binned)
    
    npix = 12*nside^2
    # now we create the interpolation
    nodes = (1:npix, log10.(En_arr))
    itp_ = Interpolations.interpolate(nodes, log10.(map_for_itp), (NoInterp(),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    itp = extrapolate(itp_, (Throw(), Line()))
    exposure_map(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    
    return exposure_map, En_arr
end


function get_exposure_map_interpolation_v3(jld2_artifact::JLD2Artifact)

    nside = jld2_artifact.nside
    filepath = joinpath(MapGenData.fits_cache, "gtexpcube2.h5")

    Emin_micro, Emax_micro = MapGenData.get_E_bins()

    

    map_binned = MapGenData._read_exposure_map_helper(nside, filepath, Emin_micro, Emax_micro, Emin_micro, Emax_micro) 
    map_for_itp = convert(Matrix{Float64}, map_binned)
    
    function exposure_map(pix::Int, En::Energy) 
        @assert Emin_micro[1] <= ustrip(u"MeV", En) <= Emax_micro[end] "Energy out of range"
        idx = searchsortedlast(Emin_micro, ustrip(u"MeV",En))
        
        return map_for_itp[pix,idx]
    end

    function exposure_map(pix::Vector{Int}, En::Energy) 
        @assert Emin_micro[1] <= ustrip(u"MeV", En) <= Emax_micro[end] "Energy out of range"
        idx = searchsortedlast(Emin_micro, ustrip(u"MeV",En))
        
        return map_for_itp[pix,idx]
    end
    
    return exposure_map
end

#%%
nside=1024
Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])
jld2_artifact = JLD2Artifact("./", nside, Emin_macro, Emax_macro)

expmap_itp = MapGenData.get_exposure_map_interpolation(jld2_artifact)

@btime expmap_itp(1,1u"GeV")

model_heal, energy_fg1 = MapGenData.read_galactic_fg_v07(jld2_artifact)

model_heal

i=1
Emin = jld2_artifact.Emin_array[i]*u"MeV"
Emax = jld2_artifact.Emax_array[i]*u"MeV"
lmax=nside*4
PSF_theta = MapGenData.get_PSF_theta(jld2_artifact)
            
PSF_theta_bini(theta) = PSF_theta(theta, i)
PSF_theta_bini(0.2432)
model_smoothed1, energy_fg1_filtered1 = MapGenData.convolve_fg_model_with_PSF(model_heal, lmax, PSF_theta_bini, energy_fg1; Emin=Emin, Emax=Emax)

# model_smoothed1 = load("gf_newversion_v7.jld2")["gf_s"]

model_smoothed1[1,1]

# model_smoothed[1,1]
# energy_fg1_filtered
# #%% counts map
# # now we interpolate the model in energy and store the interpolations in an array
gf_model_interpolated1 = Vector{Function}(undef, length(jld2_artifact.Emin_array))


gf_model_interpolated1[1], _ = MapGenData.interpolate_galactic_fg(model_smoothed1, energy_fg1_filtered1, jld2_artifact)


# # now we compute the integral of the model in each energy bin in order to obtain the foreground model in units of counts
gf_integral1 = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))

exp_map_E1 = get_exposure_map_interpolation_v3(jld2_artifact::JLD2Artifact)
#MapGenData.get_exposure_map_interpolation(jld2_artifact)
@info "Convolving the smoothed foreground template with the exposure map and computing the integral in each energy bin"
p = Progress(length(gf_integral1))
for i in eachindex(gf_integral1)
    Emin = jld2_artifact.Emin_array[i]*u"MeV"
    Emax = jld2_artifact.Emax_array[i]*u"MeV"
    gf_integral1[i] = MapGenData.process_galactic_fg_smoothed_counts(gf_model_interpolated1[i], exp_map_E1, jld2_artifact; Emin=Emin, Emax=Emax)
    next!(p)
end


gf_integral1[1][1]
#%% comparison with results
gf_oldversion_v7 = load("gf_newversion_v7.jld2")["gf"]
gf_oldversion_v7[1]
#%% using old code

nside=1024
Earr = [1000, 10_000]*u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])
jld2_artifact = JLD2Artifact("./", nside, Emin_macro, Emax_macro)
PSF_theta_ = MapGenData.get_PSF_theta(jld2_artifact)
            
PSF_theta2(theta) = PSF_theta_(theta, 1)

@memoize function W_beam_fermi(l::Int)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta2(theta)
    return min(2*pi*quadgk(arg, 0, deg2rad(29), rtol=1e-5)[1],1)
end

function apply_W_beam(map::AbstractVector, lmax::Int)
    map_ = HealpixMap{Float64, RingOrder}(map)
    return apply_W_beam(map_, lmax)
end

function apply_W_beam(map::HealpixMap{T,O}, lmax::Int) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
    W_beam_arr = zeros(lmax+1) 
    @info "Computing the beam window function"
    p = Progress(lmax+1)
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi(l)
        next!(p)
    end
    almxfl!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end

function map_coordinates(input, grid)
    (dec2,ra2) = grid
    (lenx, leny) = size(input)
    x = range(minimum(dec2), maximum(dec2), length=lenx)
    y = range(minimum(ra2), maximum(ra2), length=leny)
    knots = (x,y)
    itp = linear_interpolation(knots, input)
    return itp.(dec2, ra2)
end

function read_galactic_fg_v07(fg_folder, nside::Int)
    res = Resolution(nside)
    npix = 12 * nside * nside

    # galactic diffuse map
    file = FITS("$(fg_folder)/gll_iem_v07.fits", "r")
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
        model_heal[:,ii] .= map_coordinates(model[ii,:,:],[dec2,ra2]) 
    end

    return model_heal, energy_fg1
end

function convolve_fg_model_with_PSF(model_heal::AbstractMatrix, lmax::Int)
    model_heal_smoothed = zeros(size(model_heal))
    @showprogress "Applying PSF to the FG template" for i in 1:size(model_heal)[2]
        model_heal_smoothed[:,i] .= apply_W_beam(view(model_heal,:,i), lmax) 
    end
    return model_heal_smoothed
end

function convolve_fg_model_with_PSF(model_heal::AbstractMatrix, lmax::Int, energy_fg1::AbstractArray; Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    emin_ = ustrip(u"MeV", Emin)
    emax_ = ustrip(u"MeV", Emax)
    filter_ = (1:length(energy_fg1))[emin_.<=energy_fg1.<=emax_]
    filter = max(minimum(filter_)-1, 1):min(maximum(filter_)+1,length(energy_fg1))

    model_heal_smoothed = convolve_fg_model_with_PSF(model_heal[:,filter], lmax)
    model_heal_smoothed[model_heal_smoothed.<0] .= 1e-100
    return model_heal_smoothed, energy_fg1[filter]
end

function interpolate_galactic_fg(model_heal_smoothed::AbstractArray, energy_fg1_filtered::AbstractArray, nside::Int)
    nchannels = length(energy_fg1_filtered)
    npix = 12*nside^2

    if nside != 1024
        map_for_itp = zeros(npix, nchannels)
        for ch in 1:nchannels
            map_hp = HealpixMap{Float64, RingOrder}(model_heal_smoothed[:,ch])
            map_for_itp[:,ch] .= udgrade(map_hp, nside)
        end
    else
        map_for_itp = model_heal_smoothed
    end
    
    # now create an inteprolation in E for ease of use
    nodes = (collect(1:npix), log10.(energy_fg1_filtered))
    itp = Interpolations.interpolate(nodes, log10.(map_for_itp), (Gridded(Constant()),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    gf_model_interpolated(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    return gf_model_interpolated, energy_fg1_filtered*u"MeV"
end

function read_exposure_map_E(exppath, nside::Int, kind::String="SV")
    @assert nside <= 2048 "nside must be <= 2048"
    @assert ispow2(nside) "nside must be a power of 2"

    if kind =="SV"
        filepath = "$exppath/gtexpcube2.fits"
    # elseif kind =="UCV"
    #     filepath = "$exppath/UCV/w9w745_ULTRACLEANVETO_t1_nside2048_expcube.fits"
    else
        @error "Kind $kind not supported. Choose `SV` or `UCV`"
    end
    nchannels = 10 
    map_ = zeros(12*2048^2, nchannels)
    FITS(filepath, "r") do f
        for channel in 1:nchannels
            map_[:,channel] .= read(f[2],"ENERGY$channel")
        end
    end
    
    En_arr = uconvert.(u"MeV", 10 .^ range(log10(1), log10(10),length=10).*u"GeV")
    
    if nside == 2048
        map_for_itp = map_
        npix = nside2npix(2048)
    else
        npix = nside2npix(nside)
        map_for_itp = zeros(npix, nchannels)
        for ch in 1:nchannels
            map_hp = HealpixMap{Float64, RingOrder}(map_[:,ch])
            map_for_itp[:,ch] .= udgrade(map_hp, nside)
        end
    end
    
    nodes = (collect(1:npix), log10.(ustrip.(u"MeV", En_arr)))
    # itp = Interpolations.interpolate(nodes, log10.(map_for_itp), (Gridded(Constant()),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    itp = Interpolations.interpolate(nodes, log10.(map_for_itp), (NoInterp(),Gridded(Linear()))) 
    exposure_map(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    
    return exposure_map, En_arr
end

function process_galactic_fg_smoothed_counts(gf_model_interpolated::Function, exppath::String, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV", kind::String="SV")
    exp_map_E, _ = read_exposure_map_E(exppath, nside, kind)
    npix = 12*nside^2
    srpixel = 4pi/npix
    
    gf_integral = HealpixMap{Float64, RingOrder}(nside)
    p = Progress(npix)
    @threads for pix in 1:npix
        arg(En) = gf_model_interpolated(pix, En)*exp_map_E(pix, En)*srpixel
        gf_integral[pix] = ustrip(u"MeV", quadgk(arg, Emin, Emax, rtol=1e-4)[1])
        next!(p)
    end
    return gf_integral
end
#%%

nside=1024
model_heal, energy_fg1 = MapGenData.read_galactic_fg_v07(jld2_artifact)
model_heal

lmax=4*nside
model_heal_smoothed, energy_fg1_filtered = convolve_fg_model_with_PSF(model_heal, lmax, energy_fg1; Emin=1u"GeV", Emax=10u"GeV")
model_heal_smoothed
energy_fg1_filtered
model_heal_smoothed[1,1]

PSF_theta2(0.2432)
#%%
exppath = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_1_10_GeV_w9-745_xgam"
Emin=1u"GeV"
Emax=10u"GeV"
itp, _ = interpolate_galactic_fg(model_heal_smoothed, energy_fg1_filtered, nside)
gf_integral = process_galactic_fg_smoothed_counts(itp, exppath, nside, Emin, Emax, "SV")

# #%%
# gf_integral[1]
#%%
function pl_f(Agal::Float64, Fiso::Float64, fore_map::AbstractArray, k_arr::AbstractArray, fac::AbstractArray)

    lambda_arr = Agal*fore_map .+ Fiso.* fac # valore previsto 
    llarg = -k_arr .* log.(lambda_arr) .+ lambda_arr
    return sum(llarg[.!isnan.(llarg)])
end

function fit_foreground(;Agal_0=1.0, Fiso_0=1e-7, limit_Agal=[0.2,1.8], limit_Fiso=[1e-8,1e-6], fermi_map, GF, exposure_map, srpixel, pls_mask)
    fac = view(exposure_map.pixels .* srpixel, pls_mask)
    fore_map = view(GF.pixels, pls_mask)
    k_arr = view(fermi_map, pls_mask) # valore osservato

    _f(Agal, Fiso) = pl_f(Agal, Fiso, fore_map, k_arr, fac)

    m = Minuit(_f, Agal_0, Fiso_0)
    m.limits[1] = limit_Agal
    m.limits[2] = limit_Fiso
    m.errordef = Minuit.LIKELIHOOD
    m.migrad()
    m.hesse()

    return [(m.values[1], m.errors[1]), (m.values[2], m.errors[2])]
end
#%%
npix = 12*1024^2
srpixel=4pi/npix
data = load("data_1024.jld2")
fermi_map = data["fermi_map"][1]
pls_mask = data["pls_mask"]
exposure_map = data["exposure_map"][1]

# GF = load("gf_newversion_v7.jld2")["gf"]
GF = load("gf_oldversion_v7.jld2")["gf"]
# GF = gf_integral1[1] 
#sum(load("gf_newversion_v7.jld2")["gf_s"], dims=2)[:].*exposure_map.*srpixel
sum(GF)
GF2 = load("gf_newversion_v7.jld2")["gf"]
sum(GF2)
#%%
fit_foreground(;Agal_0=1.0, Fiso_0=1e-7, limit_Agal=[0.2,1.8], limit_Fiso=[1e-8,1e-6], fermi_map, GF, exposure_map, srpixel, pls_mask)
#%%
MapGenData.W_beam_fermi(100, PSF_theta2)
W_beam_fermi(100)

exp_map_E, En = read_exposure_map_E(exppath, 1024, "SV")

exp_map_E2, En2 = MapGenData.get_exposure_map_interpolation(jld2_artifact)

exp_map_E(3,1u"GeV")

exp_map_E2(3,1u"GeV")

En
En2
#%%
filepath = "$exppath/gtexpcube2.fits"
file = FITS(filepath, "r")

Enarra = 10 .^ range(log10(1000), log10(10000),length=11)
Emean = sqrt.(Enarra[1:end-1].*Enarra[2:end])

Enfile=read(file["ENERGIES"], "Energy")

Enfile[1]
En[1]
En2[1]

Emean[1]