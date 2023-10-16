# function W_beam_fermi(l::Int, PSF_theta::Function)
#     arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
#     return min(2*pi*quadgk(arg, 0, deg2rad(19), rtol=1e-5)[1],1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
# end

@memoize function W_beam_fermi(l::Int, PSF_theta::Function)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
    np = max(1000, l)
    i1 = quad(arg, 0, deg2rad(2), method=:gausslegendre, order=np)[1]
    i2 = quad(arg, deg2rad(2), deg2rad(19), method=:gausslegendre, order=np)[1]
    return min(2*pi*(i1+i2),1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
end

@memoize function get_W_beam_arr(lmax::Int, PSF_theta::Function)
    W_beam_arr = zeros(lmax+1) 
    p = Progress(lmax, desc="W_beam:", enabled=false)
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi(l, PSF_theta)
        next!(p)
    end
    return W_beam_arr
end

function apply_W_beam(map::AbstractVector, lmax::Int, PSF_theta::Function)
    map_ = HealpixMap{Float64, RingOrder}(map)
    return apply_W_beam(map_, lmax, PSF_theta)
end

function apply_W_beam(map::HealpixMap{T,O}, lmax::Int, PSF_theta::Function) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
    W_beam_arr = get_W_beam_arr(lmax, PSF_theta)
    almxfl!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end

function apply_W_beam(map::HealpixMap{T,O}, W_beam_arr::AbstractVector) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
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

 
function _convolve_fg_model_with_PSF_helper(model_heal::AbstractMatrix, lmax::Int, PSF_theta::Function; verbose=false)
    model_heal_smoothed = zeros(size(model_heal))
    W_beam_arr = get_W_beam_arr(lmax, PSF_theta)
    p = Progress(size(model_heal)[2], desc="PSF:", enabled=verbose)
    for i in 1:size(model_heal)[2]
        model_heal_smoothed[:,i] .= apply_W_beam(view(model_heal,:,i), W_beam_arr) 
        next!(p)
    end
    return model_heal_smoothed
end

function convolve_fg_model_with_PSF(model_heal::AbstractMatrix, lmax::Int, PSF_theta::Function, energy_fg1::AbstractArray; Emin::Energy, Emax::Energy)
    emin_ = ustrip(u"MeV", Emin)
    emax_ = ustrip(u"MeV", Emax)
    filter_ = (1:length(energy_fg1))[emin_.<=energy_fg1.<=emax_]
    filter = max(minimum(filter_)-1, 1):min(maximum(filter_)+1,length(energy_fg1))
    model_heal_filtered = view(model_heal,:,filter)
    model_heal_smoothed = _convolve_fg_model_with_PSF_helper(model_heal_filtered, lmax, PSF_theta)
    model_heal_smoothed[model_heal_smoothed.<0] .= 1e-100
    return model_heal_smoothed, energy_fg1[filter]
end

function downgrade_smoothed_template(jld2_artifact::JLD2Artifact, hres_path::String, nside::Int; verbose=true)
    dict_hr = load(hres_path)
    model_heal, energy_fg1 = dict_hr["gf_v07"], dict_hr["E"]
    # we now downgrade the maps and then export them
    model_heal_downgraded = Vector{Matrix{Float64}}(undef, length(jld2_artifact.Emin_array))
    
    steps = 0
    for i in eachindex(energy_fg1)
        steps += length(energy_fg1[i])
    end
    p = Progress(steps, desc="Downgrading:", enabled=verbose)

    @threads for i in eachindex(jld2_artifact.Emin_array)
        model_downgraded_i = zeros(12*nside^2, length(energy_fg1[i]))
        for j in eachindex(energy_fg1[i])
            model_downgraded_i[:,j] = ud_grade(HealpixMap{Float64, RingOrder}(model_heal[i][:,j]), nside).pixels
            next!(p)
        end
        model_heal_downgraded[i] = model_downgraded_i
    end
    return model_heal_downgraded, energy_fg1
end


function write_gf_v07_map_smoothed_as_jld2(jld2_artifact::JLD2Artifact; compress=true, overwrite=false, verbose=true)
    nside = jld2_artifact.nside
    @info "Processing galactic foreground v7 map at nside=$nside"
    # npix = 12*nside^2
    PSF_theta = get_PSF_theta(jld2_artifact)
    hres_path = "$(artifact_cache)/galactic_foreground_v07_nside1024.jld2"
    outpath = "$(artifact_cache)/galactic_foreground_v07_nside$(nside).jld2"

    if ispath(outpath) && !overwrite
        @info "File galactic_foreground_v07_nside$(nside).jld2 already exists, skipping"
        return outpath
    end

    if ispath(hres_path) && nside < 1024

        model_heal_downgraded, energy_fg1 = downgrade_smoothed_template(jld2_artifact, hres_path, nside)
        
        dict_ = Dict{String, Any}()
        dict_["gf_v07"] = model_heal_downgraded
        dict_["E"] = energy_fg1
        save(outpath, dict_, compress=compress)
        return outpath
    else
        if nside < 1024
            @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
            It is advisable to export the galactic foreground template at nside=1024, which will allow for downgrading of the high-res map."
        end
        model_heal, energy_fg1 = read_galactic_fg_v07(jld2_artifact)
        lmax=3*nside-1
        # here
        model_heal_smoothed = Vector{Matrix{Float64}}(undef, length(jld2_artifact.Emin_array))
        energy_fg1_filtered = Vector{Vector{Float64}}(undef, length(jld2_artifact.Emin_array))
        
        @info "Convolving the foreground template with the PSF"
        p = Progress(length(jld2_artifact.Emin_array), desc="Conv:", enabled=verbose)
        for i in eachindex(jld2_artifact.Emin_array)
            Emin = jld2_artifact.Emin_array[i]*u"MeV"
            Emax = jld2_artifact.Emax_array[i]*u"MeV"
            PSF_theta_bini(theta) = PSF_theta(theta, i)
            model_heal_smoothed[i], energy_fg1_filtered[i] = convolve_fg_model_with_PSF(model_heal, lmax, PSF_theta_bini, energy_fg1; Emin=Emin, Emax=Emax)
            next!(p)
        end
        dict_ = Dict{String, Any}()
        dict_["gf_v07"] = model_heal_smoothed
        dict_["E"] = energy_fg1_filtered
        save(outpath, dict_, compress=compress)
        return outpath
    end
end

function interpolate_galactic_fg(model::AbstractArray, energy::AbstractArray, jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    # nchannels = length(energy_fg1_filtered)
    npix = 12*nside^2

    # now create an inteprolation in E for ease of use
    nodes = (1:npix, log10.(energy))
    itp_ = Interpolations.interpolate(nodes, log10.(model), (NoInterp(),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    itp = extrapolate(itp_, (Throw(), Line()))
    gf_model_interpolated(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    return gf_model_interpolated, energy*u"MeV"
end

function process_galactic_fg_smoothed_counts(gf_model_interpolated::Function, exp_map_E::Function, jld2_artifact::JLD2Artifact; Emin::Energy, Emax::Energy, verbose=true)
    nside = jld2_artifact.nside
    npix = 12*nside^2
    srpixel = 4pi/npix
    
    gf_integral = HealpixMap{Float64, RingOrder}(nside)
    p = Progress(npix, enabled=verbose)
    @threads for pix in 1:npix
        arg(En) = gf_model_interpolated(pix, En)*exp_map_E(pix, En)*srpixel
        gf_integral[pix] = ustrip(u"MeV", quadgk(arg, Emin, Emax, rtol=1e-4)[1])
        next!(p)
    end
    return gf_integral
end

function _write_gf_v07_counts_map_as_jld2_helper(outdir::String, jld2_artifact::JLD2Artifact; compress=true, verbose=true)
    # first we load the galactic foreground model that we computed previously
    nside = jld2_artifact.nside
    @info "Processing galactic foreground v7 map at nside=$nside"
    smoothed_gf_path = "$artifact_cache/galactic_foreground_v07_nside$(nside).jld2"
    data_smoothed_gf = load(smoothed_gf_path)
    model_smoothed, energy_fg1_filtered = data_smoothed_gf["gf_v07"], data_smoothed_gf["E"]

    # now we interpolate the model in energy and store the interpolations in an array
    gf_model_interpolated = Vector{Function}(undef, length(jld2_artifact.Emin_array))

    for i in eachindex(gf_model_interpolated )
        gf_model_interpolated[i], _ = interpolate_galactic_fg(model_smoothed[i], energy_fg1_filtered[i], jld2_artifact)
    end

    # now we compute the integral of the model in each energy bin in order to obtain the foreground model in units of counts
    gf_integral = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
    
    exp_map_E, _ = get_exposure_map_interpolation(jld2_artifact)
    @info "Convolving the smoothed foreground template with the exposure map and computing the integral in each energy bin"
    p = Progress(length(gf_integral), enabled=verbose)
    for i in eachindex(gf_integral)
        Emin = jld2_artifact.Emin_array[i]*u"MeV"
        Emax = jld2_artifact.Emax_array[i]*u"MeV"
        gf_integral[i] = process_galactic_fg_smoothed_counts(gf_model_interpolated[i], exp_map_E, jld2_artifact; Emin=Emin, Emax=Emax)
        next!(p)
    end

    # now we save the foreground model
    dict_ = Dict{String, Any}()
    dict_["gf_v07"] = gf_integral
    dict_["Emin"] = jld2_artifact.Emin_array
    dict_["Emax"] = jld2_artifact.Emax_array

    outpath = joinpath(outdir, "galactic_foreground_smoothed_counts.jld2")
    save(outpath, dict_, compress=compress)
    if nside == 1024
        cp(outpath, joinpath(artifact_cache, "galactic_foreground_smoothed_counts_nside1024.jld2"), force=true)
    end
    return outpath
end

function write_gf_v07_counts_map_as_jld2(outdir::String, jld2_artifact::JLD2Artifact; compress=true)
    nside = jld2_artifact.nside
    gfpath = joinpath(artifact_cache, "galactic_foreground_smoothed_counts_nside1024.jld2")
    if nside < 1024 && isfile(gfpath)
        # we can skip the computation, and we shall downgrade the model
        data_1024 = load(gfpath)
        gf_1024 = data_1024["gf_v07"]

        data = Dict{String, Any}()
        data["Emin"] = data_1024["Emin"]
        data["Emax"] = data_1024["Emax"]

        gf_integral = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
        for i in eachindex(gf_integral)
            gf_integral[i] = ud_grade(gf_1024[i], nside, power=-2)
        end
        data["gf_v07"] = gf_integral

        outpath = joinpath(outdir, "galactic_foreground_smoothed_counts.jld2")
        save(outpath, data, compress=compress)
        return outpath
    else
        return _write_gf_v07_counts_map_as_jld2_helper(outdir, jld2_artifact; compress=compress)
    end
end

