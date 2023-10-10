using Pkg
using Revise
Pkg.activate(".")

using MapGenData
using HDF5
using Unitful
using StaticArrays
using JLD2
using Base.Threads
using Healpix
using ProgressMeter
using QuadGK
using BenchmarkTools
using MultiQuad

using PyPlot
#%%
# MapGenData.clear_cache()
#%%
@info "Using nthreads = $(nthreads())"

hdf5_folder = "/lhome/ific/a/aamerio/data/fermi/output/sourceveto_nside2048_front_0.5_1000_GeV/hdf5"
artifacts_folder = "/lhome/ific/a/aamerio/data/artifacts"
artifact_cache = MapGenData.artifact_cache

Earr = [500, 1000, 2000, 5000, 10_000, 50_000, 200_000, 1_000_000] * u"MeV"

Emin_macro = ustrip.(u"MeV", Earr[1:end-1])
Emax_macro = ustrip.(u"MeV", Earr[2:end])

#%%
function almxfl_multi!(alm::Alm{Complex{T}}, fl::AA) where {T<:Number,AA<:AbstractArray{T,1}}

    lmax = alm.lmax
    mmax = alm.mmax
    fl_size = length(fl)

    if lmax + 1 > fl_size
        fl = [fl; zeros(lmax + 1 - fl_size)]
    end
    for m = 0:mmax
        @threads for l = m:lmax
            i = almIndex(alm, l, m)
            alm.alm[i] = alm.alm[i] * fl[l+1]
        end
    end
end
#%%

#%%
nside = 1024
jld2_artifact = JLD2Artifact(artifacts_folder, nside, Emin_macro, Emax_macro)
#%%
model_heal, energy_fg1 = MapGenData.read_galactic_fg_v07(jld2_artifact)
lmax = nside * 5 ÷ 2
# here
model_heal_smoothed = Vector{Matrix{Float64}}(undef, length(jld2_artifact.Emin_array))
energy_fg1_filtered = Vector{Vector{Float64}}(undef, length(jld2_artifact.Emin_array))
#%%
i = 1
PSF_theta = MapGenData.get_PSF_theta(jld2_artifact)
Emin = jld2_artifact.Emin_array[i] * u"MeV"
Emax = jld2_artifact.Emax_array[i] * u"MeV"
PSF_theta_bini(theta) = PSF_theta(theta, i)
PSF_theta_bini(0.3)
map_= view(model_heal,:,i)
#%%
using LegendrePolynomials
function W_beam_fermi_v2(l::Int, PSF_theta::Function)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
    np = max(1000, l)
    i1 = quad(arg, 0, deg2rad(2), method=:gausslegendre, order=np)[1]
    i2 = quad(arg, deg2rad(2), deg2rad(19), method=:gausslegendre, order=np)[1]
    return min(2*pi*(i1+i2),1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
end

# model_heal_smoothed[i], energy_fg1_filtered[i] = MapGenData.convolve_fg_model_with_PSF(model_heal, lmax, PSF_theta_bini, energy_fg1; Emin=Emin, Emax=Emax)
#%%
function almxfl_multi!(alm::Alm{Complex{T}}, fl::AA) where {T<:Number,AA<:AbstractArray{T,1}}

    lmax = alm.lmax
    mmax = alm.mmax
    fl_size = length(fl)

    if lmax + 1 > fl_size
        fl = [fl; zeros(lmax + 1 - fl_size)]
    end
    for m = 0:mmax
        @threads for l = m:lmax
            i = almIndex(alm, l, m)
            alm.alm[i] = alm.alm[i] * fl[l+1]
        end
    end
end

function apply_W_beam_v2(map::HealpixMap{T,O}, lmax::Int, PSF_theta::Function) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
    W_beam_arr = zeros(lmax+1) 
    p = Progress(lmax, desc="W_beam:")
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi_v2(l, PSF_theta)
        next!(p)
    end
    println("computing multiplication")
    almxfl_multi!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end

function apply_W_beam_v2(map::AbstractVector, lmax::Int, PSF_theta::Function)
    map_ = HealpixMap{Float64, RingOrder}(map)
    return apply_W_beam_v2(map_, lmax, PSF_theta)
end
#%%
map1 = @time MapGenData.apply_W_beam(map_, 3*nside, PSF_theta_bini)
map2 = @time apply_W_beam_v2(map_, 3*nside, PSF_theta_bini)
#%%
abs_diff = abs.(map1 .- map2)./map1
println(maximum(abs_diff))
#%%
# @btime W_beam_fermi($lmax, $PSF_theta_bini)
# #%%
# l=2000
# @btime MapGenData.W_beam_fermi(l, PSF_theta_bini)
# @btime W_beam_fermi_v2(l, PSF_theta_bini)
# arg(theta) = 2*pi*sin(theta)*Pl(cos(theta), l)*PSF_theta(theta, 1)
# @btime begin
#     a1 = quad(arg, 0, deg2rad(2), method=:gausslegendre, order=2_000)[1] 
#     a2 = quad(arg, deg2rad(2),deg2rad(19), method=:gausslegendre, order=1_000)[1] 
#     a1+a2
# end

# quadgk(arg, 0, deg2rad(2), deg2rad(19), rtol=1e-3, order=20)[1]
# #%%
# theta = range(0,19,length=10000)
# y = arg.(deg2rad.(theta))
# xrad = deg2rad.(theta)
# #%%
# using PyCall
# np = pyimport("numpy")
# np.trapz(y, xrad)

# #%%
# plt.clf()
# plt.plot(theta, y)
# # plt.yscale("log")
# plt.show()
# plt.gcf()