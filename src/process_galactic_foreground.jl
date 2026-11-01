@memoize function W_beam_fermi(l::Int, PSF_theta::Function)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
    # return min(2*pi*quadgk(arg, 0, deg2rad(2), deg2rad(29), rtol=1e-5)[1],1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
    return min(2*pi*quadgk(arg, 0,  deg2rad(29), rtol=1e-5)[1],1) 
end

# @memoize function W_beam_fermi(l::Int, PSF_theta::Function)
#     arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
#     np = max(1000, l)
#     i1 = quad(arg, 0, deg2rad(2), method=:gausslegendre, order=np)[1]
#     i2 = quad(arg, deg2rad(2), deg2rad(19), method=:gausslegendre, order=np)[1]
#     return min(2*pi*(i1+i2),1) # TODO: 19 gradi potrebbero essere pochi, magari meglio suare 29
# end

@memoize function get_W_beam_arr(lmax::Int, PSF_theta::Function)
    W_beam_arr = zeros(lmax+1) 
    p = Progress(lmax, desc="W_beam:", enabled=false)
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi(l, PSF_theta)
        next!(p)
    end
    return W_beam_arr
end

function apply_W_beam(map::AbstractVector{T}, lmax::Int, PSF_theta::Function) where {T<:Real}
    map_ = HealpixMap{T, RingOrder}(map)
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
    lmax = length(W_beam_arr)-1
    alm_map = map2alm(map, lmax=lmax)
    almxfl!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end

function apply_W_beam(map::AbstractVector{T}, W_beam_arr::AbstractVector) where {T<:Real}
    map_ = HealpixMap{T, RingOrder}(map)
    return apply_W_beam(map_, W_beam_arr)
end

##### 
"""
    map_coordinates(input, grid)

Bilinear samples of `input` (`nlat` by `nlon`) at the 1-based fractional indices
`grid = (row, col)`: grid node `(i, j)` sits at index `(i, j)`. Longitude (columns) is
periodic: `col` is wrapped into `[0.5, nlon + 0.5)` and the grid is padded with one
column on each side (the opposite edge column), so a point beside the seam interpolates
between the first and last columns. Latitude (rows) is clamped to the edge rows.
"""
function map_coordinates(input, grid)
    (dec2, ra2) = grid
    (lenx, leny) = size(input)
    padded = hcat(input[:, leny], input, input[:, 1])
    itp = linear_interpolation((1:lenx, 1:leny+2), padded)
    row = clamp.(dec2, 1, lenx)
    col = mod.(ra2 .- 0.5, leny) .+ 1.5 # wrapped to [0.5, nlon + 0.5), then +1 for the pad
    return itp.(row, col)
end

"""
    car_pixel_indices(header, nside)

1-based fractional `(row, col)` of the RING pixel centres of `nside` on a galactic CAR
grid, from the `CRPIX`, `CRVAL` and `CDELT` keywords of `header` (FITS WCS: pixel `p`
of axis `k` sits at `CRVAL_k + (p - CRPIX_k) * CDELT_k`), longitude wrapped into
`(-180, 180]` around `CRVAL1`.
"""
function car_pixel_indices(header::FITSHeader, nside::Int)
    res = Resolution(nside)
    npix = nside2npix(nside)
    crpix1, crval1, cdelt1 = Float64(header["CRPIX1"]), Float64(header["CRVAL1"]), Float64(header["CDELT1"])
    crpix2, crval2, cdelt2 = Float64(header["CRPIX2"]), Float64(header["CRVAL2"]), Float64(header["CDELT2"])
    row = zeros(npix)
    col = zeros(npix)
    for j in 1:npix
        (theta, phi) = pix2angRing(res, j)
        b = 90 - rad2deg(theta)
        lon_wrapped = 180 - mod(180 - (rad2deg(phi) - crval1), 360) # in (-180, 180]
        row[j] = crpix2 + (b - crval2) / cdelt2
        col[j] = crpix1 + lon_wrapped / cdelt1
    end
    return row, col
end

"""
    sample_galactic_fg_v07(filepath, nside)

Sample every energy plane of a galactic CAR template cube (gll_iem_v07 layout: primary
image, `ENERGIES` extension with column `energy` in MeV) at the RING pixel centres of
`nside`, using the WCS keywords of the primary header. Returns `(model_heal (npix, nE),
energy_fg1)`.
"""
function sample_galactic_fg_v07(filepath::AbstractString, nside::Int)
    file = FITS(filepath, "r")
    energy_fg1 = read(file[2], "energy")
    header = read_header(file[1])
    model_ = Float64.(read(file[1]))
    close(file)

    model = permutedims(model_, (3,2,1)) # (nE, nlat, nlon)
    row, col = car_pixel_indices(header, nside)

    model_heal = ones((length(row), length(energy_fg1)))
    for ii in eachindex(energy_fg1)
        model_heal[:,ii] .= map_coordinates(model[ii,:,:], [row, col])
    end
    return model_heal, energy_fg1
end

function read_galactic_fg_v07(jld2_artifact::JLD2Artifact)
    # galactic diffuse map
    filepath = joinpath(fits_cache, "gll_iem_v07.fits")
    return sample_galactic_fg_v07(filepath, jld2_artifact.nside)
end

function read_galactic_fg_v05(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    res = Resolution(nside)
    npix = 12 * nside^2

    # galactic diffuse map
    filepath = joinpath(fits_cache, "gll_iem_v05_rev1.fit")
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
    nres = 8 # 8
    dec2 = ( -dec2 * 180 / pi .+ 180 ) .* nres
    ra2 = ( ( +ra2 * 180 / pi .+ 360 .+ 180 ) .% 360 ) .* nres #fg v5 uses a different coord convention for the grid
    # 0-based grid indices above; map_coordinates takes 1-based indices and wraps the longitude
    dec2 .+= 1
    ra2 .+= 1

    model_heal = ones((npix,eneb_fg))

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
    if length(filter_) == 0
        # we extrapolate the foreground model
        @warn "Extrapolating the foreground model to Emin=$Emin and Emax=$Emax"
        npix = size(model_heal)[1]
        knots = (1:npix, log10.(energy_fg1))
        itp_ = Interpolations.interpolate(knots, log10.(model_heal), (NoInterp(),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
        itp = extrapolate(itp_, (Throw(), Line()))
        logEn =  range(log10(emin_), log10(emax_), length=3)
        model_heal_filtered = 10 .^ itp.(reshape(collect(1:npix),:,1), reshape(collect(logEn),1,:))

        model_heal_smoothed = _convolve_fg_model_with_PSF_helper(model_heal_filtered, lmax, PSF_theta)
        model_heal_smoothed[model_heal_smoothed.<0] .= 1e-100
        return model_heal_smoothed, 10 .^logEn
    
    else
        filter = max(minimum(filter_)-1, 1):min(maximum(filter_)+1,length(energy_fg1))
        model_heal_filtered = view(model_heal,:,filter)
        model_heal_smoothed = _convolve_fg_model_with_PSF_helper(model_heal_filtered, lmax, PSF_theta)
        model_heal_smoothed[model_heal_smoothed.<0] .= 1e-100
        return model_heal_smoothed, energy_fg1[filter]
    end
end

function downgrade_smoothed_template(jld2_artifact::JLD2Artifact, hres_path::String, nside::Int; verbose=true)
    dict_hr = load(hres_path)
    model_heal, energy_fg1 = dict_hr["gf"], dict_hr["E"]
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


function write_gf_map_smoothed_as_jld2(jld2_artifact::JLD2Artifact; version=7, compress=false, overwrite=false, verbose=true)
    nside = jld2_artifact.nside
    @info "Smoothing galactic foreground v$(version) map at nside=$nside"
    # npix = 12*nside^2
    PSF_theta = get_PSF_theta(jld2_artifact)
    hres_path = "$(artifact_cache)/galactic_foreground_v0$(version)_nside1024.jld2"
    outpath = "$(artifact_cache)/galactic_foreground_v0$(version)_nside$(nside).jld2"

    if ispath(outpath) && !overwrite
        @info "File galactic_foreground_v0$(version)_nside$(nside).jld2 already exists, skipping"
        return outpath
    end

    if ispath(hres_path) && nside < 1024

        model_heal_downgraded, energy_fg1 = downgrade_smoothed_template(jld2_artifact, hres_path, nside)
        
        dict_ = Dict{String, Any}()
        dict_["gf"] = model_heal_downgraded
        dict_["E"] = energy_fg1
        save(outpath, dict_, compress=compress)
        return outpath
    else
        if nside < 1024
            @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
            It is advisable to export the galactic foreground template at nside=1024, which will allow for downgrading of the high-res map."
        end
        if version == 7
            model_heal, energy_fg1 = read_galactic_fg_v07(jld2_artifact)
        elseif version == 5
            model_heal, energy_fg1 = read_galactic_fg_v05(jld2_artifact)
        else
            error("Version $version not supported")
        end
        lmax=4*nside
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
        dict_["gf"] = model_heal_smoothed
        dict_["E"] = energy_fg1_filtered
        save(outpath, dict_, compress=compress)
        return outpath
    end
end

function interpolate_galactic_fg(model::AbstractArray, energy::AbstractArray, jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    # nchannels = length(energy_fg1_filtered)
    npix = 12*nside^2

    size(model, 1) == npix || error("model has $(size(model, 1)) pixels, expected $npix for nside $nside")
    size(model, 2) == length(energy) || error("model has $(size(model, 2)) planes for $(length(energy)) energies")
    length(energy) >= 2 || error("need at least two template planes, got $(length(energy))")

    # linear in log10 E - log10 f between the planes, and along the end segments beyond the
    # first and last plane. Done by hand: Interpolations.jl's (NoInterp(), Gridded(Linear()))
    # with (Throw(), Line()) holds the value flat outside the knots (Interpolations 0.15.1).
    logE = log10.(collect(Float64, energy))
    all(diff(logE) .> 0) || error("template energies must be strictly increasing, got $energy")
    logm = log10.(model)
    gf_model_interpolated(pix::Int, En::Energy) = 10^_loglog_line(logm, logE, pix, log10(ustrip(u"MeV", En))) * u"MeV^-1" # the gf is in units of cm^-2 s^-1 sr^-1 MeV^-1. While cm^-2 s^-1 sr^-1 disappear when multiplying by the exposure and srpix, we need to divide by MeV to get the correct units for the Energy integration
    return gf_model_interpolated, energy*u"MeV"
end

# log10 f of pixel `pix` at log10 E = x: the straight line through the two planes that
# bracket x, or through the first or last two planes when x is outside them.
function _loglog_line(logm::AbstractMatrix, logE::AbstractVector, pix::Int, x::Real)
    i = clamp(searchsortedlast(logE, x), 1, length(logE) - 1)
    y0 = logm[pix, i]
    y1 = logm[pix, i + 1]
    return y0 + (y1 - y0) * (x - logE[i]) / (logE[i + 1] - logE[i])
end

function process_galactic_fg_smoothed_counts(gf_model_interpolated::Function, exp_map_E::Function, jld2_artifact::JLD2Artifact; Emin::Energy, Emax::Energy, verbose=true)
    nside = jld2_artifact.nside
    npix = 12*nside^2
    srpixel = 4pi/npix
    
    gf_integral = HealpixMap{Float64, RingOrder}(nside)
    # the gf model is in units of cm^-2 s^-1 sr^-1 MeV^-1, the exposure is in units of cm^2 s, srpixel is the conversion factor between sr and pixel. As a result, the integral is in units of counts per pixel
    p = Progress(npix, enabled=verbose)
    @threads for pix in 1:npix
        arg(En) = gf_model_interpolated(pix, En)*exp_map_E(pix, En)*srpixel
        gf_integral[pix] = quadgk(arg, Emin, Emax, rtol=1e-4)[1]
        next!(p)
    end
    return gf_integral
end

function _write_gf_counts_map_as_jld2_helper(outdir::String, jld2_artifact::JLD2Artifact; version=7, compress=false, verbose=true)
    # first we load the galactic foreground model that we computed previously
    nside = jld2_artifact.nside
    @info "Convolving the galactic foreground map at nside=$nside with the exposure map"
    smoothed_gf_path = "$artifact_cache/galactic_foreground_v0$(version)_nside$(nside).jld2"
    data_smoothed_gf = load(smoothed_gf_path)
    model_smoothed, energy_fg1_filtered = data_smoothed_gf["gf"], data_smoothed_gf["E"]

    # now we interpolate the model in energy and store the interpolations in an array
    gf_model_interpolated = Vector{Function}(undef, length(jld2_artifact.Emin_array))

    for i in eachindex(gf_model_interpolated )
        gf_model_interpolated[i], _ = interpolate_galactic_fg(model_smoothed[i], energy_fg1_filtered[i], jld2_artifact)
    end

    # now we compute the integral of the model in each energy bin in order to obtain the foreground model in units of counts
    gf_integral = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
    
    exp_map_E = get_exposure_map_interpolation(jld2_artifact)
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
    dict_["gf"] = gf_integral
    dict_["Emin"] = jld2_artifact.Emin_array
    dict_["Emax"] = jld2_artifact.Emax_array

    outpath = joinpath(outdir, "galactic_foreground_v0$(version)_smoothed_counts.jld2")
    save(outpath, dict_, compress=compress)
    if nside == 1024
        cp(outpath, joinpath(artifact_cache, "galactic_foreground_v0$(version)_smoothed_counts_nside1024.jld2"), force=true)
    end
    return outpath
end

function write_gf_counts_map_as_jld2(outdir::String, jld2_artifact::JLD2Artifact; version=7, compress=false)
    nside = jld2_artifact.nside
    gfpath = joinpath(artifact_cache, "galactic_foreground_v0$(version)_smoothed_counts_nside1024.jld2")
    if nside < 1024 && isfile(gfpath)
        @info "Downgrading the galactic foreground counts map from nside=1024 to nside=$nside"
        # we can skip the computation, and we shall downgrade the model
        data_1024 = load(gfpath)
        gf_1024 = data_1024["gf"]

        data = Dict{String, Any}()
        data["Emin"] = data_1024["Emin"]
        data["Emax"] = data_1024["Emax"]

        gf_integral = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
        for i in eachindex(gf_integral)
            gf_integral[i] = ud_grade(gf_1024[i], nside, power=-2)
        end
        data["gf"] = gf_integral

        outpath = joinpath(outdir, "galactic_foreground_v0$(version)_smoothed_counts.jld2")
        save(outpath, data, compress=compress)
        return outpath
    else
        return _write_gf_counts_map_as_jld2_helper(outdir, jld2_artifact; version=version, compress=compress)
    end
end

