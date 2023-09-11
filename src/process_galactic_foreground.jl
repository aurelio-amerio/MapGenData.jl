# const fg_folder = "$(fermitools_path)/galactic_foreground"
# const fg_v7_path = "$fg_folder/gll_iem_v07.fits"
# const fg_v5_path = "$fg_folder/gll_iem_v05.fits"

function W_beam_fermi(l::Int)
    arg(theta) = sin(theta)*Pl(cos(theta), l)*PSF_theta(theta)
    return min(2*pi*quadgk(arg, 0, deg2rad(29), rtol=1e-5)[1],1)
end

function apply_W_beam(map::AbstractVector, lmax::Int)
    map_ = HealpixMap{Float64, RingOrder}(map)
    return apply_W_beam(map_, lmax)
end

function apply_W_beam(map::HealpixMap{T,O}, lmax::Int) where {T<:AbstractFloat, O<:RingOrder}
    alm_map = map2alm(map, lmax=lmax)
    W_beam_arr = zeros(lmax+1) 
    @threads for l in 0:lmax
        W_beam_arr[l+1] = W_beam_fermi(l)
    end
    almxfl!(alm_map, W_beam_arr)
    map_smoothed = alm2map(alm_map, npix2nside(length(map)))
    return map_smoothed
end


function download_gll_iem_v07(fg_folder)
    @info "Downloading gll_iem_v07"

    outpath = "$(fg_folder)/gll_iem_v07.fits"
    url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits"

    p_ = Progress(1)

    function progress_bar(total, current)
        p_.n=total
        update!(p_, current)
    end

    download(url, outpath,progress=progress_bar)
end

function download_gll_iem_v05(fg_folder)
    @info "Downloading gll_iem_v05"
    outpath = "$(fg_folder)/gll_iem_v05.fits"
    url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v05_rev1.fit"

    p_ = Progress(1)

    function progress_bar(total, current)
        p_.n=total
        update!(p_, current)
    end

    download(url, outpath,progress=progress_bar)
end

function download_smoothed_gf_v7(dir)
    @info "Downloading smoothed galactic foreground v7"
    id = "1BbYELQrqYIPuHP3S51RX-dmBphgq5g9Z"
    file = "galactic_foreground_v07_nside1024_1-10GeV.zip"

    outpath = "$dir/$file"
    download_from_gdrive(id, outpath)
end


# function check_foregrounds()
#     @info "Checking galactic foreground files"
#     if !(isfile(fg_v7_path))
#         download_gll_iem_v07()
#     end

#     if !(isfile(fg_v5_path))
#         download_gll_iem_v05()
#     end
# end

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

function read_galactic_fg_v05(fg_folder, nside::Int)
    res = Resolution(nside)
    npix = 12 * nside * nside

    # galactic diffuse map
    file = FITS("$(fg_folder)/gll_iem_v05.fits", "r")
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
    ra2 = ( ( -ra2 * 180 / pi .+ 360 .+ 180 ) .% 360 ) .* nres

    model_heal = ones((npix,eneb_fg))

    ra2[ra2.>2879].=2879
    dec2[dec2.>1440].=1440

    for ii in eachindex(energy_fg1)
        model_heal[:,ii] .= map_coordinates(model[ii,:,:],[dec2,ra2]) 
    end

    return model_heal, energy_fg1
end

# function process_galactic_fg_v07_old(fg_folder, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    
#     model_heal, energy_fg1 = read_galactic_fg_v07(fg_folder, nside::Int)

#     Emin = ustrip(u"MeV", Emin) 
#     Emax = ustrip(u"MeV", Emax) 

#     map_gf = HealpixMap{Float64, RingOrder}(nside)
#     p = Progress(npix)
#     @threads for j in 1:npix
#         interfunc_tmp = linear_interpolation(log10.(energy_fg1), log10.(model_heal[j,:]))

#         interfunc(E) = 10 ^ interfunc_tmp(log10(E))

#         map_gf[j] = quadgk(interfunc, Emin, Emax, rtol=1e-3)[1]
#         next!(p)
#     end
#     return map_gf
# end

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

function write_gf_v07_map_smoothed_as_jld2(fg_folder, dir, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    @info "Processing galactic foreground v7 map"
    # processing the foreground is expensive. 
    if nside < 512
        @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
        It is advisable to export the galactic foreground template at nside=1024 first, to allow for map downgrading."
    end
    model_heal, energy_fg1 = read_galactic_fg_v07(fg_folder, nside)
    lmax=4*nside
    model_heal_smoothed, energy_fg1_filtered = convolve_fg_model_with_PSF(model_heal, lmax, energy_fg1; Emin=Emin, Emax=Emax)
    dict_ = Dict{String, Array}()
    dict_["gf_v07"] = model_heal_smoothed
    dict_["E"] = energy_fg1_filtered
    save("$dir/galactic_foreground_v07_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV.jld2", dict_,compress=true)
    return
end

function write_gf_v05_map_smoothed_as_jld2(fg_folder, dir, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
    @info "Processing galactic foreground v5 map"
    # processing the foreground is expensive. 
    if nside < 512
        @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
        It is advisable to export the galactic foreground template at nside=1024 first, to allow for map downgrading."
    end
    model_heal, energy_fg1 = read_galactic_fg_v05(fg_folder, nside)
    lmax=4*nside
    model_heal_smoothed, energy_fg1_filtered = convolve_fg_model_with_PSF(model_heal, lmax, energy_fg1; Emin=Emin, Emax=Emax)
    dict_ = Dict{String, Array}()
    dict_["gf_v05"] = model_heal_smoothed
    dict_["E"] = energy_fg1_filtered
    # If we have already computed the HR fg at nside 1024, we will just downgrade it
    # if nside<1024 && ispath(path1024)
    #     data = load(path1024)
    #     map_1024 = HealpixMap{Float64, RingOrder}(data["gf_v07"])
    #     dict_ = Dict{String, Array}()
    #     dict_["gf_v07"] = MapGen.ud_grade(map_1024, nside, power=-2)
    # else
    #     if nside < 512
    #         @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
    #         It is advisable to export the galactic foreground template at nside=1024 first, to allow for map downgrading."
    #     end
    #     model_heal, energy_fg1 = read_galactic_fg_v07(fg_folder, nside)
    #     model_heal_smoothed, energy_fg1_filtered = convolve_fg_model_with_PSF(model_heal, lmax, energy_fg1; Emin=Emin, Emax=Emax)
    #     dict_ = Dict{String, Array}()
    #     dict_["gf_v07"] = model_heal_smoothed
    #     dict_["E"] = energy_fg1_filtered
    # end
    save("$dir/galactic_foreground_v05_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV.jld2", dict_,compress=true)
    return
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

# #todo edit
function write_gf_v07_map_as_jld2(exppath, dir, basepath, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV", kind::String="SV")
    @info "Processing galactic foreground v7 map"
    path1024 = "$basepath/galactic_foreground_smoothed_v7/galactic_foreground_v07_nside1024_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV.jld2"
    # processing the foreground is expensive. 
    # If we have already computed the HR fg at nside 1024, we will just downgrade it
    
    data_1024 = load(path1024)
    itp, _ = interpolate_galactic_fg(data_1024["gf_v07"], data_1024["E"], nside)
    gf_integral = process_galactic_fg_smoothed_counts(itp, exppath, nside, Emin, Emax, kind)
    dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
    dict_["gf_v07"] = gf_integral
    
    save("$dir/gf_smoothed_counts_v07_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV_$(kind).jld2", dict_,compress=true)
    return
end


# function write_gf_v05_map_as_jld2(fg_folder, exppath, dir, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV", kind::String="SV")
#     @info "Processing galactic foreground v5 $kind map"
#     path1024 = "$dir/galactic_foreground_v05_nside1024_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV_$(kind).jld2"
#     # processing the foreground is expensive. 
#     # If we have already computed the HR fg at nside 1024, we will just downgrade it
#     if nside<1024 && ispath(path1024)
#         map_1024 = load(path1024)["gf_v05"]
#         dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
#         dict_["gf_v05"] = MapGen.ud_grade(map_1024, nside, power=-2)
#     else
#         if nside < 512
#             @warn "Computing the smoothed galactic foreground directly at nside=$nside, results may be inaccurate.
#             It is advisable to export the galactic foreground template at nside=1024 first, to allow for map downgrading."
#         end
#         model_heal, energy_fg1 = read_galactic_fg_v05(fg_folder, nside)
#         gf_model_interpolated, energy_fg1 = interpolate_galactic_fg(model_heal, energy_fg1, nside; Emin=Emin, Emax=Emax, lmax=4*nside)
#         map_ = process_galactic_fg_smoothed_counts(gf_model_interpolated, exppath, nside, Emin, Emax, kind)
#         dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
#         dict_["gf_v05"] = map_
#     end
#     save("$dir/galactic_foreground_v05_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV_$(kind).jld2", dict_,compress=compress)
#     return
# end

# TODO

# function process_galactic_fg_v05_old(fg_folder, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
#     res = Resolution(nside)
#     npix = 12 * nside * nside

#     # galactic diffuse map
#     file = FITS("$(fg_folder)/gll_iem_v05.fits", "r")

#     energy_fg1 = read(file[2], "Energy")
#     model_ = Float64.(read(file[1]))
#     close(file)

#     model = permutedims(model_, (3,2,1))
#     eneb_fg = length(energy_fg1)

#     dec2 = zeros(npix)
#     ra2 = zeros(npix)
#     for j in 1:npix
#         (dec2[j], ra2[j]) = pix2angRing(res, j)
#     end
#     nres = 8 # 8
#     dec2 = ( -dec2 * 180 / pi .+ 180 ) .* nres
#     ra2 = ( ( +ra2 * 180 / pi .+ 360 .+ 180 ) .% 360 ) .* nres

#     model_heal = ones((npix,eneb_fg))

#     ra2[ra2.>2879].=2879
#     dec2[dec2.>1440].=1440

#     for ii in eachindex(energy_fg1)
#         model_heal[:,ii] .= map_coordinates(model[ii,:,:],[dec2,ra2]) 
#     end

#     Emin = ustrip(u"MeV", Emin) 
#     Emax = ustrip(u"MeV", Emax) 

#     map_gf = HealpixMap{Float64, RingOrder}(nside)
#     p = Progress(npix)
#     @threads for j in 1:npix
#         interfunc_tmp = linear_interpolation(log10.(energy_fg1), log10.(model_heal[j,:]))

#         interfunc(E) = 10 ^ interfunc_tmp(log10(E))

#         map_gf[j] = quadgk(interfunc, Emin, Emax, rtol=1e-3)[1]
#         next!(p)
#     end
#     return map_gf
# end

# function write_gf_v05_map_as_jld2(fg_folder, dir, nside::Int, Emin::Energy=1u"GeV", Emax::Energy=10u"GeV")
#     @info "Processing galactic foreground v5 map"
#     map_ = process_galactic_fg_v05(fg_folder, nside, Emin, Emax)
#     dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
#     dict_["gf_v05"] = map_
#     save("$dir/galactic_foreground_v05_nside$(nside)_$(ustrip(u"GeV",Emin))-$(ustrip(u"GeV",Emax))GeV.jld2", dict_, compress=compress)
#     return
# end