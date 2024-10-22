function make_artifacts(fits_artifact::FITSArtifact, jld2_artifact::JLD2Artifact; skip_fits=false, skip_jld2=false)
    fetch_fermilat_data(fits_artifact)
    if ! skip_fits
        make_fits_artifact(fits_artifact)
    end
    if ! skip_jld2
        make_jld2_artifacts(jld2_artifact)
    end
    return
end

function fetch_fermilat_data(fits_artifact::FITSArtifact)
    @info "Copying .h5 fits files"
    #first we copy all the files to the scratch directory, 
    #then we copy them to the temp directory to make the artifact
    if ! isfile(joinpath(fits_cache, "gtbin.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtbin.h5"), joinpath(fits_cache, "gtbin.h5"), force=true)
    end
    if ! isfile(joinpath(fits_cache, "gtexpcube2.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtexpcube2.h5"), joinpath(fits_cache, "gtexpcube2.h5"), force=true)
    end
    if ! isfile(joinpath(fits_cache, "gtpsf.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtpsf.h5"), joinpath(fits_cache, "gtpsf.h5"), force=true)
    end

    @info "Fetching gll_psc"
    #4FGL-DR4 catalog
    gll_psc_url_32 = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/gll_psc_v32.fit"

    #4FGL-DR3 catalog
    gll_psc_url_30 = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/12yr_catalog/gll_psc_v30.fit"
    gll_psc_url_31 = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/12yr_catalog/gll_psc_v31.fit"

    if ! isfile(joinpath(fits_cache, "gll_psc_v30.fit"))
        Downloads.download(gll_psc_url_30, joinpath(fits_cache,"gll_psc_v30.fit"))
    end

    if ! isfile(joinpath(fits_cache, "gll_psc_v31.fit"))
        Downloads.download(gll_psc_url_31, joinpath(fits_cache,"gll_psc_v31.fit"))
    end
    
    if ! isfile(joinpath(fits_cache, "gll_psc_v32.fit"))
        Downloads.download(gll_psc_url_32, joinpath(fits_cache,"gll_psc_v32.fit"))
    end

    @info "Fetching foreground template v7"
    #foreground template
    if ! isfile(joinpath(fits_cache, "gll_iem_v07.fits"))
        if isfile(joinpath(fits_artifact.fitsdir, "gll_iem_v07.fits"))
            cp(joinpath(fits_artifact.fitsdir, "gll_iem_v07.fits"), joinpath(fits_cache, "gll_iem_v07.fits"), force=true)
        else
            gll_iem_url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits"
            Downloads.download(gll_iem_url, joinpath(fits_cache, "gll_iem_v07.fits"))
        end
    end

    @info "Fetching foreground template v5"
    if ! isfile(joinpath(fits_cache, "gll_iem_v05_rev1.fit"))
        if isfile(joinpath(fits_artifact.fitsdir, "gll_iem_v05_rev1.fit"))
            cp(joinpath(fits_artifact.fitsdir, "gll_iem_v05_rev1.fit"), joinpath(fits_cache, "gll_iem_v05_rev1.fit"), force=true)
        else
            gll_iem_url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v05_rev1.fit"
            Downloads.download(gll_iem_url, joinpath(fits_cache, "gll_iem_v05_rev1.fit"))
        end
    end
    return
end

# make the artifact for the fits files
function make_fits_artifact(fits_artifact::FITSArtifact)
    tmp_dir = mktempdir()

    # copy from the artifact cache to the tmp directory
    # cp(joinpath(fits_cache, "gtbin.h5"), joinpath(tmp_dir, "gtbin.h5"), force=true)
    # cp(joinpath(fits_cache, "gtexpcube2.h5"), joinpath(tmp_dir, "gtexpcube2.h5"), force=true)
    # cp(joinpath(fits_cache, "gtpsf.h5"), joinpath(tmp_dir, "gtpsf.h5"), force=true)

    cp(joinpath(fits_cache, "gll_psc_v30.fit"), joinpath(tmp_dir, "gll_psc_v30.fit"), force=true)
    cp(joinpath(fits_cache, "gll_psc_v31.fit"), joinpath(tmp_dir, "gll_psc_v31.fit"), force=true)
    cp(joinpath(fits_cache, "gll_psc_v32.fit"), joinpath(tmp_dir, "gll_psc_v32.fit"), force=true)

    # cp(joinpath(fits_cache, "gll_iem_v07.fits"), joinpath(tmp_dir, "gll_iem_v07.fits"), force=true)
    # cp(joinpath(fits_cache, "gll_iem_v05_rev1.fit"), joinpath(tmp_dir, "gll_iem_v05_rev1.fit"), force=true)

    mkpath(fits_artifact.outdir)

    @info "Creating tarball"
    fits_artifact_id = artifact_from_directory(tmp_dir)
    sha256_fits = Pkg.archive_artifact(fits_artifact_id, "$(fits_artifact.outdir)/fits.tar.gz")
 
    return "$(fits_artifact.outdir)/fits.tar.gz", "SHA256: $sha256_fits"
end

function make_jld2_artifacts(jld2_artifact::JLD2Artifact)
    tmp_dir = mktempdir()
    tmp_dir_fg_v5 = mktempdir()
    tmp_dir_fg_v7 = mktempdir()
    nside = jld2_artifact.nside
    cache = mkpath(joinpath(artifact_cache, "nside$nside"))

    @info "Processing jld2 artifacts for nside $nside"
    @info "Exporting Fermi map"
    if ! isfile(joinpath(cache,"fermi_map.jld2"))
        write_fermi_map_as_jld2(cache, jld2_artifact; compress=true)
    end
    cp(joinpath(cache,"fermi_map.jld2"), joinpath(tmp_dir, "fermi_map.jld2"), force=true)

    @info "Exporting exposure map"
    if ! isfile(joinpath(cache,"exposure_map.jld2"))
        write_exposure_map_as_jld2(cache, jld2_artifact; compress=true)
    end
    cp(joinpath(cache,"exposure_map.jld2"), joinpath(tmp_dir, "exposure_map.jld2"), force=true)

    @info "Exporting PSF"
    if ! isfile(joinpath(cache,"PSF.jld2"))
        write_PSF_as_jld2(cache, jld2_artifact; compress=true)
    end
    cp(joinpath(cache,"PSF.jld2"), joinpath(tmp_dir, "PSF.jld2"), force=true)

    @info "Exporting galactic foreground"
    # make sure we have already computed the foreground at nside 1024
    # if ! isfile(joinpath(artifact_cache, "galactic_foreground_v7_nside1024.jld2"))
    #     let
    #         @info "Computing smoothed galactic foreground v07 at nside 1024"
    #         jld2_artifact_tmp = JLD2Artifact("./", 1024, jld2_artifact.Emin_array, jld2_artifact.Emax_array)
    #         write_gf_map_smoothed_as_jld2(jld2_artifact_tmp, version=7)
    #     end
    # end

    # if ! isfile(joinpath(artifact_cache, "galactic_foreground_v05_nside1024.jld2"))
    #     let
    #         @info "Computing smoothed galactic foreground v5 at nside 1024"
    #         jld2_artifact_tmp = JLD2Artifact("./", 1024, jld2_artifact.Emin_array, jld2_artifact.Emax_array)
    #         write_gf_map_smoothed_as_jld2(jld2_artifact_tmp, version=5)
    #     end
    # end

    # now we compute the gf as a counts map
    if ! isfile(joinpath(cache,"galactic_foreground_v07_smoothed_counts.jld2"))
        write_gf_map_smoothed_as_jld2(jld2_artifact, version=7)
        write_gf_counts_map_as_jld2(cache, jld2_artifact; version=7, compress=false)
    end
    cp(joinpath(cache,"galactic_foreground_v07_smoothed_counts.jld2"), joinpath(tmp_dir_fg_v7, "galactic_foreground_v07_smoothed_counts.jld2"), force=true)

    if ! isfile(joinpath(cache,"galactic_foreground_v05_smoothed_counts.jld2"))
        write_gf_map_smoothed_as_jld2(jld2_artifact, version=5)
        write_gf_counts_map_as_jld2(cache, jld2_artifact; version=5, compress=false)
    end
    cp(joinpath(cache,"galactic_foreground_v05_smoothed_counts.jld2"), joinpath(tmp_dir_fg_v5, "galactic_foreground_v05_smoothed_counts.jld2"), force=true)

    @info "Creating tarball"
    jld2_artifact_id = artifact_from_directory(tmp_dir)
    sha256_jld2 = Pkg.archive_artifact(jld2_artifact_id, "$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz")

    fg7_artifact_id = artifact_from_directory(tmp_dir_fg_v7)
    sha256_fg7 = Pkg.archive_artifact(fg7_artifact_id, "$(jld2_artifact.outdir)/jld2_data_fg_v7_n$nside.tar.gz")

    fg5_artifact_id = artifact_from_directory(tmp_dir_fg_v5)
    sha256_fg5 = Pkg.archive_artifact(fg5_artifact_id, "$(jld2_artifact.outdir)/jld2_data_fg_v5_n$nside.tar.gz")

    return #"$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz", "SHA256: $sha256_jld2"
end

function clear_cache(;clear_fermilat_data=false)
    delete_scratch!(MapGenData, artifact_cache_name)
    if clear_fermilat_data
        delete_scratch!(MapGenData, "fits_cache")
    end
    @info "Cache cleared"
    __init__() # reinitialize the cache
end
