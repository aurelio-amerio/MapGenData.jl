function W_beam_fermi(l::Int, PSF_theta::Function)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
    return min(2*pi*quadgk(arg, 0, deg2rad(19), rtol=1e-5)[1],1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
end

function apply_W_beam(map::AbstractVector, lmax::Int, PSF_theta::Function)
    map_ = HealpixMap{Float64, RingOrder}(map)
    return apply_W_beam(map_, lmax, PSF_theta)
end

function apply_W_beam(map::HealpixMap{T,O}, lmax::Int, PSF_theta::Function) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
    W_beam_arr = zeros(lmax+1) 
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi(l, PSF_theta)
    end
    almxfl!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end

##### 
function map_coordinates(input, grid)
    (dec2,ra2) = grid
    (lenx, leny) = size(input)
    x = range(minimum(dec2), maximum(dec2), length=lenx)
    y = range(minimum(ra2), maximum(ra2), length=leny)
    knots = (x,y)
    itp = linear_interpolation(knots, input)
    return itp.(dec2, ra2)
end

function read_galactic_fg_v07(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    res = Resolution(nside)
    npix = 12 * nside^2

    # galactic diffuse map
    filepath = joinpath(artifact_cache, "fits", "gll_iem_v07.fits")

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
        model_heal[:,ii] .= map_coordinates(model[ii,:,:],[dec2,ra2]) 
    end

    return model_heal, energy_fg1
end

 
function _convolve_fg_model_with_PSF_helper(model_heal::AbstractMatrix, lmax::Int, PSF_theta::Function)
    model_heal_smoothed = zeros(size(model_heal))
    @showprogress "Applying PSF to the FG template" for i in 1:size(model_heal)[2]
        model_heal_smoothed[:,i] .= apply_W_beam(view(model_heal,:,i), lmax, PSF_theta) 
    end
    return model_heal_smoothed
end

function convolve_fg_model_with_PSF(model_heal::AbstractMatrix, lmax::Int, PSF_theta::Function, energy_fg1::AbstractArray; Emin::Energy, Emax::Energy)
    emin_ = ustrip(u"MeV", Emin)
    emax_ = ustrip(u"MeV", Emax)
    filter_ = (1:length(energy_fg1))[emin_.<=energy_fg1.<=emax_]
    filter = max(minimum(filter_)-1, 1):min(maximum(filter_)+1,length(energy_fg1))

    model_heal_smoothed = _convolve_fg_model_with_PSF_helper(model_heal[:,filter], lmax, PSF_theta)
    model_heal_smoothed[model_heal_smoothed.<0] .= 1e-100
    return model_heal_smoothed, energy_fg1[filter]
end


function write_gf_v07_map_smoothed_as_jld2(jld2_artifact::JLD2Artifact, outdir::String; compress=true)
    @info "Processing galactic foreground v7 map"
    nside = jld2_artifact.nside
    PSF_theta = get_PSF_theta(jld2_artifact)
    
    # processing the foreground is expensive. 
    if nside < 512
        @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
        It is advisable to export the galactic foreground template at nside=1024 first, to allow for map downgrading."
    end
    model_heal, energy_fg1 = read_galactic_fg_v07(jld2_artifact)
    lmax=4*nside
    # here
    model_heal_smoothed = Vector{Matrix{Float64}}(undef, length(jld2_artifact.Emin_array))
    energy_fg1_filtered = Vector{Vector{Float64}}(undef, length(jld2_artifact.Emin_array))
    for i in eachindex(jld2_artifact.Emin_array)
        Emin = jld2_artifact.Emin_array[i]*u"MeV"
        Emax = jld2_artifact.Emax_array[i]*u"MeV"
        PSF_theta_bini(theta) = PSF_theta(theta, i)
        model_heal_smoothed[i], energy_fg1_filtered[i] = convolve_fg_model_with_PSF(model_heal, lmax, PSF_theta_bini, energy_fg1; Emin=Emin, Emax=Emax)
    end
    dict_ = Dict{String, Any}()
    dict_["gf_v07"] = model_heal_smoothed
    dict_["E"] = energy_fg1_filtered
    save("$outdir/galactic_foreground_v07_nside$(nside).jld2", dict_, compress=compress)
    return "$outdir/galactic_foreground_v07_nside$(nside).jld2"
end
#=
function interpolate_galactic_fg(model::AbstractArray, energy::AbstractArray, jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    # nchannels = length(energy_fg1_filtered)
    npix = 12*nside^2

    # if nside != 1024
    #     map_for_itp = zeros(npix, nchannels)
    #     for ch in 1:nchannels
    #         map_hp = HealpixMap{Float64, RingOrder}(model_heal_smoothed[:,ch])
    #         map_for_itp[:,ch] .= udgrade(map_hp, nside)
    #     end
    # else
    #     map_for_itp = model_heal_smoothed
    # end
    
    # now create an inteprolation in E for ease of use
    nodes = (collect(1:npix), log10.(energy))
    itp = Interpolations.interpolate(nodes, log10.(model), (Gridded(Constant()),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    gf_model_interpolated(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    return gf_model_interpolated, energy*u"MeV"
end

function process_galactic_fg_smoothed_counts(gf_model_interpolated::Function, jld2_artifact::JLD2Artifact; Emin::Energy, Emax::Energy)
    nside = jld2_artifact.nside
    exp_map_E, _ = get_exposure_map_interpolation(jld2_artifact)
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

# #todo edit
function write_gf_v07_map_as_jld2(jld2_artifact::JLD2Artifact, outdir::String; compress=true)
    @info "Processing galactic foreground v7 map"
    model, energy = read_galactic_fg_v07(jld2_artifact)
    
    itp, _ = interpolate_galactic_fg(model_heal_smoothed::AbstractArray, energy_fg1_filtered::AbstractArray, jld2_artifact::JLD2Artifact)
    gf_integral = process_galactic_fg_smoothed_counts(itp, exppath, nside, Emin, Emax, kind)
    dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
    dict_["gf_v07"] = gf_integral
    
    save("$dir/gf_smoothed_counts_v07_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV_$(kind).jld2", dict_,compress=true)
    return
end
=#