# make the artifact for the fits files
function make_fits_artifact(fits_artifact::FITSArtifact)
    tmp_dir = mktempdir()
    mkpath(joinpath(artifact_cache, "fits"))
    # tmp_dir = mkpath("/tmp/tmp_art")
    @info "Copying .h5 fits files"
    #first we copy all the files to the scratch directory, 
    #then we copy them to the temp directory to make the artifact
    if ! isfile(joinpath(artifact_cache, "gtbin.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtbin.h5"), joinpath(artifact_cache, "fits", "gtbin.h5"), force=true)
    end
    if ! isfile(joinpath(artifact_cache, "gtexpcube2.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtexpcube2.h5"), joinpath(artifact_cache, "fits", "gtexpcube2.h5"), force=true)
    end
    if ! isfile(joinpath(artifact_cache, "gtpsf.h5"))
        cp(joinpath(fits_artifact.fitsdir, "gtpsf.h5"), joinpath(artifact_cache, "fits", "gtpsf.h5"), force=true)
    end

    # copy from the artifact cache to the tmp directory
    cp(joinpath(artifact_cache, "fits", "gtbin.h5"), joinpath(tmp_dir, "gtbin.h5"), force=true)
    cp(joinpath(artifact_cache, "fits", "gtexpcube2.h5"), joinpath(tmp_dir, "gtexpcube2.h5"), force=true)
    cp(joinpath(artifact_cache, "fits", "gtpsf.h5"), joinpath(tmp_dir, "gtpsf.h5"), force=true)
    

    @info "Fetching gll_psc"
    #4FGL-DR4 catalog
    gll_psc_url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/gll_psc_v32.fit"
    
    if ! isfile(joinpath(artifact_cache, "fits", "gll_psc_v32.fit"))
        Downloads.download(gll_psc_url, joinpath(artifact_cache, "fits","gll_psc_v32.fit"))
    end

    cp(joinpath(artifact_cache, "fits", "gll_psc_v32.fit"), joinpath(tmp_dir, "gll_psc_v32.fit"), force=true)

    @info "Fetching foreground template"
    #foreground template
    if ! isfile(joinpath(artifact_cache, "fits", "gll_iem_v07.fits"))
        if isfile(joinpath(fits_artifact.fitsdir, "gll_iem_v07.fits"))
            cp(joinpath(fits_artifact.fitsdir, "gll_iem_v07.fits"), joinpath(artifact_cache, "fits", "gll_iem_v07.fits"), force=true)
        else
            gll_iem_url = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits"
            Downloads.download(gll_iem_url, joinpath(artifact_cache, "fits", "gll_iem_v07.fits"))
        end
    end

    cp(joinpath(artifact_cache, "fits", "gll_iem_v07.fits"), joinpath(tmp_dir, "gll_iem_v07.fits"), force=true)

    mkpath(fits_artifact.outdir)

    @info "Creating tarball"
    fits_artifact_id = artifact_from_directory(tmp_dir)
    sha256_fits = Pkg.archive_artifact(fits_artifact_id, "$(fits_artifact.outdir)/fits.tar.gz")
 
    return "$(fits_artifact.outdir)/fits.tar.gz", "SHA256: $sha256_fits"
end

function make_jld2_artifacts(jld2_artifact::JLD2Artifact)
    tmp_dir = mktempdir()
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
    if ! isfile(joinpath(artifact_cache, "galactic_foreground_v07_nside$(nside).jld2"))
        jld2_artifact_tmp = JLD2Artifact(artifacts_cache, nside, jld2_artifact.Emin_array, jld2_artifact.Emax_array)
        write_gf_v07_map_smoothed_as_jld2(jld2_artifact_tmp)
    end

    if ! isfile(joinpath(cache,"galactic_foreground_smoothed_counts.jld2"))
        write_gf_v07_map_smoothed_as_jld2(jld2_artifact)
        write_gf_v07_counts_map_as_jld2(cache, jld2_artifact; compress=true)
    end
    cp(joinpath(cache,"galactic_foreground_smoothed_counts.jld2"), joinpath(tmp_dir, "galactic_foreground_smoothed_counts.jld2"), force=true)

    @info "Creating tarball"
    jld2_artifact_id = artifact_from_directory(tmp_dir)
    sha256_jld2 = Pkg.archive_artifact(jld2_artifact_id, "$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz")
    return "$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz", "SHA256: $sha256_jld2"
end

function clear_cache()
    delete_scratch!(MapGenData, "artifact_cache")
    @info "Cache cleared"
end
