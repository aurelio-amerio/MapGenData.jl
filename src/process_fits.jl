function check_bins(Emin_micro::AbstractVector, Emax_micro::AbstractVector, Emin_macro::AbstractVector, Emax_macro::AbstractVector)
    flag = true
    @assert all(Emin_micro .< Emax_micro) "Emin_micro must be < Emax_micro"
    @assert all(Emin_macro .< Emax_macro) "Emin_macro must be < Emax_macro"

    flag = flag & all(Emax_macro .<= maximum(Emax_micro)) # out of bounds

    flag = flag & all(Emin_macro .>= minimum(Emin_micro)) # out of bounds

    flag = flag & all([Emin in Emin_micro for Emin in Emin_macro])
    flag = flag & all([Emax in Emax_micro for Emax in Emax_macro])
    return flag
end

function get_E_bins()
    filepath = joinpath(artifact_cache, "fits", "gtbin.h5")
    @assert isfile(filepath) "gtbin.h5 not found. Run `make_fits_artifact` first"
    fid = h5open(filepath, "r") 
    Emin_ = fid["EBOUNDS/DATA/E_MIN"][:] #keV
    Emax_ = fid["EBOUNDS/DATA/E_MAX"][:] #keV
    close(fid)
    return Emin_ / 1000, Emax_ / 1000 #MeV
end

function get_E_bins(units::EnergyFreeUnits)
    Emin_, Emax_ = get_E_bins() .* u"MeV"
    return uconvert.(units, Emin_), uconvert.(units, Emax_)
end


function read_fermi_map(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    filepath = joinpath(artifact_cache, "fits", "gtbin.h5")
    @assert isfile(filepath) "gtbin.h5 not found. Run `make_fits_artifact` first"

    Emin_micro, Emax_micro = get_E_bins()

    fid = h5open(filepath, "r") 

    # initialize the map
    map_binned = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
    for i in eachindex(map_binned)
        map_binned[i] = HealpixMap{Float64, RingOrder}(2048)
    end
    
    #fill the map
    @showprogress "Reading map" for i in eachindex(jld2_artifact.Emin_array)
        for j in eachindex(Emin_micro)
            if (Emin_micro[j] >= jld2_artifact.Emin_array[i]) & (Emax_micro[j] <= jld2_artifact.Emax_array[i])
                map_binned[i].pixels .+= fid["SKYMAP/DATA"]["CHANNEL$j"][:]
            end
        end
    end
    close(fid)

    if nside == 2048
        return map_binned
    else
        map_binned_ud = Vector{HealpixMap{Float64, RingOrder}}(undef, length(jld2_artifact.Emin_array))
        @showprogress "Downgrading map" for i in eachindex(map_binned_ud)
            map_binned_ud[i] = ud_grade(map_binned[i], nside, power=-2)
        end
        return map_binned_ud
    end
end

function write_fermi_map_as_jld2(outdir, jld2_artifact::JLD2Artifact; compress=true)
    @info "Processing Fermi photon counts map"
    map_binned = read_fermi_map(jld2_artifact)
    dict_ = Dict{String, Any}()
    dict_["fermi_map"] = map_binned
    dict_["Emin"] = jld2_artifact.Emin_array # MeV
    dict_["Emax"] = jld2_artifact.Emax_array # MeV

    save("$outdir/fermi_map.jld2", dict_, compress=compress)
    return "$outdir/fermi_map.jld2"
end

"""
    _read_exposure_map_helper(filepath::String, Emin_micro::AbstractArray, Emax_micro::AbstractArray, Emin_macro::AbstractArray, Emax_macro::AbstractArray)

Reads the exposure map from the gtexpcube2.fits file and returns a vector of Healpix maps, one for each energy bin.

# Arguments
- `filepath`: path to the gtexpcube2.fits file
- `Emin_micro`: array of lower energy bin edges in MeV, as contained in the fits file
- `Emax_micro`: array of upper energy bin edges in MeV, as contained in the fits file
- `Emin_macro`: array of lower energy bin edges in MeV, as per analysys requirements
- `Emax_macro`: array of upper energy bin edges in MeV, as per analysys requirements
"""
function _read_exposure_map_helper(nside::Int, filepath::String, Emin_micro::AbstractArray, Emax_micro::AbstractArray, Emin_macro::AbstractArray, Emax_macro::AbstractArray)
    fid = h5open(filepath, "r") 

    # initialize the map
    map_binned = Vector{HealpixMap{Float64, RingOrder}}(undef, length(Emin_macro))
    for i in eachindex(map_binned)
        map_binned[i] = HealpixMap{Float64, RingOrder}(nside)
    end
    
    #fill the map
    bin_counts = zeros(Int, length(Emin_macro))
    map_cache = HealpixMap{Float64, RingOrder}(2048)
    @showprogress "Reading map" for i in eachindex(Emin_macro)
        for j in eachindex(Emin_micro)
            if (Emin_micro[j] >= Emin_macro[i]) & (Emax_micro[j] <= Emax_macro[i])
                map_cache.pixels .= fid["HPXEXPOSURES/DATA"]["ENERGY$j"][:]
                map_binned[i].pixels .+= ud_grade(map_cache, nside)
                bin_counts[i] += 1
            end
        end
    end
    close(fid)

    for i in eachindex(map_binned)
        map_binned[i].pixels ./= bin_counts[i]
    end

    return map_binned
end

function read_exposure_map(jld2_artifact::JLD2Artifact)
    nside = jld2_artifact.nside
    filepath = joinpath(artifact_cache, "fits", "gtexpcube2.h5")

    Emin_micro, Emax_micro = get_E_bins()

    map_binned = _read_exposure_map_helper(nside, filepath, Emin_micro, Emax_micro, jld2_artifact.Emin_array, jld2_artifact.Emax_array)

    return map_binned
end

function write_exposure_map_as_jld2(outdir, jld2_artifact::JLD2Artifact; compress=true)
    @info "Processing exposure map"
    map_binned = read_exposure_map(jld2_artifact)
    dict_ = Dict{String, Any}()
    dict_["exposure_map"] = map_binned
    dict_["Emin"] = jld2_artifact.Emin_array # MeV
    dict_["Emax"] = jld2_artifact.Emax_array # MeV

    save("$outdir/exposure_map.jld2", dict_, compress=compress)
    return "$outdir/exposure_map.jld2"
end

function get_PSF_arrays(jld2_artifact::JLD2Artifact)
    filepath = joinpath(artifact_cache, "fits", "gtpsf.h5")

    fid = h5open(filepath, "r")

    En=fid["PSF/DATA"]["Energy"][:] # MeV
    PSF_mat=fid["PSF/DATA"]["Psf"][:,:]
    theta=fid["THETA/DATA/Theta"][:] # deg

    close(fid)

    function PSF_mean(theta_i)
        num = zeros(length(jld2_artifact.Emin_array))
        den = zeros(length(jld2_artifact.Emin_array))
        for k in eachindex(jld2_artifact.Emin_array)
            for j in eachindex(En)
                if jld2_artifact.Emin_array[k] <= En[j] < jld2_artifact.Emax_array[k]
                    Ej = En[j] ^ (-2.4)
                    num[k] += PSF_mat[theta_i, j] * Ej
                    den[k] += Ej
                end
            end
        end
        return num ./ den
    end
    PSF_matrix = zeros(length(theta), length(jld2_artifact.Emax_array))
    for i in eachindex(theta)
        PSF_matrix[i,:] .= PSF_mean(i)
    end

    return theta, PSF_matrix
end

function get_PSF_theta(jld2_artifact::JLD2Artifact)
    theta, PSF_matrix = get_PSF_arrays(jld2_artifact)
    nbins = size(PSF_matrix)[2]
    if nbins > 1
        bin_arr = 1:nbins
        nodes = (theta, bin_arr)
        itp_ = Interpolations.interpolate(nodes, log10.(PSF_matrix), (Gridded(Linear()),NoInterp()))
        PSF_theta(theta, bin::Int) = 10 .^ itp_(rad2deg(theta), bin)
        return PSF_theta
    else
        nodes = (theta,)
        itp_ = Interpolations.interpolate(nodes, log10.(vec(PSF_matrix)), Gridded(Linear()))
        function PSF_theta_1d(theta, bin::Int) 
            @assert bin==1 "bin must be 1 if only one energy bin is present"
            return 10 .^ itp_(rad2deg(theta))
        end
        return PSF_theta_1d
    end
end

function write_PSF_as_jld2(outdir, jld2_artifact::JLD2Artifact; compress=true)
    @info "Processing PSF" 
    theta, PSF = get_PSF_arrays(jld2_artifact)
    dict_ = Dict{String, Any}()
    dict_["theta"] = theta
    dict_["PSF"] = PSF
    dict_["Emin"] = jld2_artifact.Emin_array # MeV
    dict_["Emax"] = jld2_artifact.Emax_array # MeV
    save("$outdir/PSF.jld2", dict_, compress=compress)
    return "$outdir/PSF.jld2"
end

function get_exposure_map_interpolation(jld2_artifact::JLD2Artifact)

    nside = jld2_artifact.nside
    filepath = joinpath(artifact_cache, "fits", "gtexpcube2.h5")

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




