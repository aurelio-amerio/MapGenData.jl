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
        jld2_artifact_tmp = JLD2Artifact(artifacts_folder, nside, "tmp", jld2_artifact.Emin_array, jld2_artifact.Emax_array)
        write_gf_v07_map_smoothed_as_jld2(jld2_artifact_tmp)
    end

    if ! isfile(joinpath(cache,"galactic_foreground_smoothed_counts.jld2"))
        write_gf_v07_map_smoothed_as_jld2(jld2_artifact)
        write_gf_v07_counts_map_as_jld2(artifact_cache,jld2_artifact; compress=true)
    end
    cp(joinpath(cache,"galactic_foreground_smoothed_counts.jld2"), joinpath(tmp_dir, "galactic_foreground_smoothed_counts.jld2"), force=true)

    @info "Creating tarball"
    jld2_artifact_id = artifact_from_directory(tmp_dir)
    sha256_jld2 = Pkg.archive_artifact(jld2_artifact_id, "$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz")
    return "$(jld2_artifact.outdir)/jld2_data_n$nside.tar.gz", "SHA256: $sha256_jld2"
end


function delete_cache()
    delete_scratch!(MapGenData, "artifact_cache")
end

#=

# dest = mkpath("artifacts")
# This script is used to create the artifacts for the MapGen.jl package

# fits_tmp = "tmp_art" #mktempdir()
# MapGen.download_fits(fits_tmp)

# fits_artifact_id = artifact_from_directory(fits_tmp)
# sha256_fits = Pkg.archive_artifact(fits_artifact_id, "artifacts/fits.tar.gz")
# println("artifacts/fits.tar.gz", "SHA256: $sha256_fits")

# import the PSF
# function get_PSF(years::Int=12, weeks::Int=745)
#     filepath = "$(fits_path)/$(years)y/psf/w9w$(weeks)_SV_t1_PSF.jld2"
#     data_dict = load(filepath)
#     PSF_array = data_dict["PSF"]
#     theta_array = data_dict["theta"]
#     itp = linear_interpolation(theta_array, PSF_array)

#     function PSF_itp(theta::Real)
#         return itp(theta * 180 / pi)
#     end
#     return PSF_itp
# end

# PSF_theta = get_PSF()

# now we create the artifacts for everything up to nside=1024, to modify once we have the artifacts created

# jld2_hash_n1024 = artifact_hash("jld2_data_n1024", artifact_toml)
# jld2_hash_n512 = artifact_hash("jld2_data_n512", artifact_toml)
# jld2_hash_n256 = artifact_hash("jld2_data_n256", artifact_toml)
# jld2_hash_n128 = artifact_hash("jld2_data_n128", artifact_toml)
# jld2_hash_n64 = artifact_hash("jld2_data_n64", artifact_toml)

function make_jld2_artifacts(base_dir::String, out_dir::String, name::String, nside::Int)
    fg_folder = "$base_dir/12y/galactic_foreground"
    binpath = "$base_dir/12y/output_gtbin"
    exppath = "$base_dir/12y/output_gtexpcube2"

    artifact_id = mktempdir() do artifact_dir
        @info "Processing jld2 artifacts for nside $nside"

        dir = mkpath("$artifact_dir/12y")
        MapGen.write_fermi_map_as_jld2(binpath, dir, nside, "SV")
        MapGen.write_exposure_map_as_jld2(exppath, dir, nside, "SV")
        MapGen.write_fermi_map_as_jld2(binpath, dir, nside, "UCV")
        MapGen.write_exposure_map_as_jld2(exppath, dir, nside, "UCV")
        MapGen.write_gf_v07_map_as_jld2(exppath, dir, joinpath(base_dir, "12y"), nside, 1u"GeV", 10u"GeV", "SV")
        MapGen.write_gf_v07_map_as_jld2(exppath, dir, joinpath(base_dir, "12y"), nside, 1u"GeV", 10u"GeV", "UCV")

        # MapGen.write_gf_v05_map_as_jld2(fg_folder, exppath, dir, nside, 1u"GeV", 10u"GeV", "SV")
        # MapGen.write_gf_v05_map_as_jld2(fg_folder, exppath, dir, nside, 1u"GeV", 10u"GeV", "UCV")
        artifact_from_directory(artifact_dir)
    end
    destfile = "$out_dir/$name.tar.gz"
    @info "Archiving artifact $name to $destfile"
    sha256 = Pkg.archive_artifact(artifact_id, destfile)
    return destfile, sha256
end

# function upload_jld2_artifacts_to_gist(base_dir::String, nside::Int)
#     fg_folder = "$base_dir/12y/galactic_foreground"
#     binpath = "$base_dir/12y/output_gtbin"
#     exppath = "$base_dir/12y/output_gtexpcube2"

#     artifact_id = mktempdir() do artifact_dir
#         @info "Processing jld2 artifacts for nside $nside"

#         dir = mkpath("$artifact_dir/12y")
#         MapGen.write_fermi_map_as_jld2(binpath, dir, nside, "SV")
#         MapGen.write_exposure_map_as_jld2(exppath, dir, nside, "SV")
#         MapGen.write_fermi_map_as_jld2(binpath, dir, nside, "UCV")
#         MapGen.write_exposure_map_as_jld2(exppath, dir, nside, "UCV")
#         MapGen.write_gf_v07_map_as_jld2(exppath, dir, joinpath(base_dir, "12y"), nside, 1u"GeV", 10u"GeV", "SV")
#         MapGen.write_gf_v07_map_as_jld2(exppath, dir, joinpath(base_dir, "12y"), nside, 1u"GeV", 10u"GeV", "UCV")

#         # MapGen.write_gf_v05_map_as_jld2(fg_folder, exppath, dir, nside, 1u"GeV", 10u"GeV", "SV")
#         # MapGen.write_gf_v05_map_as_jld2(fg_folder, exppath, dir, nside, 1u"GeV", 10u"GeV", "UCV")
#         artifact_from_directory(artifact_dir)
#     end
#     gist = upload_to_gist(artifact_id)
#     return gist
# end

# make_jld2_artifacts(fits_tmp, "artifacts", "jld2_data_n1024", 1024)
# make_jld2_artifacts(fits_tmp, "artifacts", "jld2_data_n512", 512)
# make_jld2_artifacts(fits_tmp, "artifacts", "jld2_data_n256", 256)
# make_jld2_artifacts(fits_tmp, "artifacts", "jld2_data_n128", 128)
# make_jld2_artifacts(fits_tmp, "artifacts", "jld2_data_n64", 64)

# gist_1024 = upload_jld2_artifacts_to_gist(fits_tmp, 1024)

# add_artifact!("Artifacts.toml", "jld2_data_n1024", gist_1024)

#%%
add_artifact!(
    "Artifacts.toml",
    "fits",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/fits.tar.gz",
    force=true,
)

add_artifact!(
    "Artifacts.toml",
    "jld2_data_n1024",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/jld2_data_n1024.tar.gz",
    force=true,
)
add_artifact!(
    "Artifacts.toml",
    "jld2_data_n512",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/jld2_data_n512.tar.gz",
    force=true,
)
add_artifact!(
    "Artifacts.toml",
    "jld2_data_n256",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/jld2_data_n256.tar.gz",
    force=true,
)
add_artifact!(
    "Artifacts.toml",
    "jld2_data_n128",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/jld2_data_n128.tar.gz",
    force=true,
)
add_artifact!(
    "Artifacts.toml",
    "jld2_data_n64",
    "https://github.com/aurelio-amerio/MapGen-data/releases/download/v0.1-beta/jld2_data_n64.tar.gz",
    force=true,
)

using Pkg
using Artifacts
Pkg.ensure_artifact_installed("jld2_data_n1024","Artifacts.toml")
artifact"jld2_data_n1024"

=#