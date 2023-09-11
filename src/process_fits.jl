function read_fermi_map(binpath, nside::Int, kind::String="SV")
    @assert nside <= 2048 "nside must be <= 2048"
    @assert ispow2(nside) "nside must be a power of 2"
    if kind == "SV"
        filepath = "$binpath/SV/w9w765_SOURCEVETO_t1_nside2048_outofbin.fits"
    elseif  kind == "UCV"
        filepath = "$binpath/UCV/w9w765_ULTRACLEANVETO_t1_nside2048_outofbin.fits"
    else
        @error "kind $kind not supported. Choose `SV` or `UCV`"
    end

    map_ = HealpixMap{Float64, RingOrder}(2048)

    FITS(filepath, "r") do f
        for channel in 1:10
            map_.pixels .+= read(f[2],"CHANNEL$channel")
        end
    end
    if nside == 2048
        return map_
    else
        return ud_grade(map_, nside, power=-2)
    end
end

function write_fermi_map_as_jld2(binpath, dir, nside::Int, kind="SV")
    @info "Processing Fermi photon counts map"
    map_ = read_fermi_map(binpath, nside, kind)
    dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
    dict_["fermi_map"] = map_
    save("$dir/w9w765_$(kind)_t1_nside$(nside)_fermi_map.jld2", dict_, compress=compress)
    return
end

function read_exposure_map(exppath, nside::Int, kind::String="SV")
    @assert nside <= 2048 "nside must be <= 2048"
    @assert ispow2(nside) "nside must be a power of 2"

    if kind =="SV"
        filepath = "$exppath/SV/w9w765_SOURCEVETO_t1_nside2048_expcube.fits"
    elseif kind =="UCV"
        filepath = "$exppath/UCV/w9w765_ULTRACLEANVETO_t1_nside2048_expcube.fits"
    else
        @error "Kind $kind not supported. Choose `SV` or `UCV`"
    end

    map_ = HealpixMap{Float64, RingOrder}(2048)

    FITS(filepath, "r") do f
        for channel in 1:10
            map_.pixels .+= read(f[2],"ENERGY$channel")
        end
    end
    map_.pixels ./= 10
    if nside == 2048
        return map_
    else
        return MapGen.ud_grade(map_, nside)
    end
end

function read_exposure_map_E(exppath, nside::Int, kind::String="SV")
    @assert nside <= 2048 "nside must be <= 2048"
    @assert ispow2(nside) "nside must be a power of 2"

    if kind =="SV"
        filepath = "$exppath/SV/w9w765_SOURCEVETO_t1_nside2048_expcube.fits"
    elseif kind =="UCV"
        filepath = "$exppath/UCV/w9w765_ULTRACLEANVETO_t1_nside2048_expcube.fits"
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
    itp = Interpolations.interpolate(nodes, log10.(map_for_itp), (Gridded(Constant()),Gridded(Linear()))) # pixel have to be exact, we interpolate linearly in energy
    exposure_map(pix::Int, En::Energy) = 10 .^ itp(pix, log10(ustrip(u"MeV", En)))
    
    return exposure_map, En_arr
end

function write_exposure_map_as_jld2(exppath, dir, nside::Int, kind="SV")
    @info "Processing exposure map"
    map_ = read_exposure_map(exppath,nside, kind)
    dict_ = Dict{String, HealpixMap{Float64, RingOrder}}()
    dict_["exposure_map"] = map_
    save("$dir/w9w765_$(kind)_t1_nside$(nside)_exposure_map.jld2", dict_, compress=compress)
    return
end

#TODO update with UCV/SV
function get_PSF_arrays(psfpath)
    filepath = "$psfpath/w9w765_SOURCEVETO_t1_psf.fits"
    f = FITS(filepath, "r") 
    En = read(f[2], "Energy")
    theta = read(f[3], "theta")
    PSF_mat = read(f[2], "Psf")
    close(f)

    function PSF_mean(theta_i)
        num = 0.
        den = 0.
        for j in eachindex(En)
            Ej = En[j] ^ (-2.4)
            num += PSF_mat[theta_i, j] * Ej
            den += Ej
        end
        return num / den
    end
    PSF_array = [PSF_mean(i) for i in eachindex(theta)]
    return theta, PSF_array
end

function write_PSF_as_jld2(psfpath, dir)
    @info "Processing PSF" 
    theta, PSF = get_PSF_arrays(psfpath)
    dict_ = Dict{String, Vector}()
    dict_["theta"] = theta
    dict_["PSF"] = PSF
    save("$dir/w9w765_SV_t1_PSF.jld2", dict_, compress=compress)
    return
end


# function get_pls_mask(dir, nside::Int)
#     pls_file_path = joinpath(dir, "mask_gll_psc_v30.jld2")
#     mask_tmp = HealpixMap{Float64, RingOrder}(Float64.(load(pls_file_path)["mask"]))

#     mask_float = udgrade(mask_tmp, nside)
#     return [pix <= 1e-12 ? false : true for pix in mask_float.pixels]
# end



