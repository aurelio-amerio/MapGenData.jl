# Standalone check of the CAR template sampler (read_galactic_fg_v07 path).
# Run: julia +1.10 --project=<MapGenData.jl> test/test_template_sampling.jl
# Builds a 64 by 33 synthetic galactic CAR cube with the same header conventions as
# gll_iem_v07 (and as MapGen-jax tests/test_preprocessing_foreground.py), samples it at
# nside 8 and compares with the fields evaluated at the pixel centres.
using MapGenData
using FITSIO
using Healpix

nlon, nlat = 64, 33
keys_ = ["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2"]
vals_ = Any["GLON-CAR", "GLAT-CAR", nlon / 2 + 0.5, (nlat + 1) / 2, 0.0, 0.0, -360.0 / nlon, 180.0 / (nlat - 1)]
header = FITSHeader(keys_, vals_, fill("", length(keys_)))

# grid node coordinates, FITS WCS: pixel p sits at CRVAL + (p - CRPIX) * CDELT
lon_nodes = [0.0 + (i - vals_[3]) * vals_[7] for i in 1:nlon]
lat_nodes = [0.0 + (j - vals_[4]) * vals_[8] for j in 1:nlat]
cube = zeros(nlon, nlat, 3) # FITS axis order: (NAXIS1 = lon, NAXIS2 = lat, NAXIS3 = energy)
for i in 1:nlon, j in 1:nlat
    cube[i, j, 1] = lat_nodes[j]                  # linear in b: bilinear is exact
    cube[i, j, 2] = cos(deg2rad(lon_nodes[i]))    # smooth in l
    cube[i, j, 3] = sin(deg2rad(lon_nodes[i]))    # changes sign across l = 180: needs the seam wrap
end

out, energies = mktempdir() do dir
    path = joinpath(dir, "synthetic_iem.fits")
    FITS(path, "w") do f
        write(f, cube; header=header)
        write(f, ["energy"], [[50.0, 100.0, 200.0]]; name="ENERGIES")
    end
    MapGenData.sample_galactic_fg_v07(path, 8)
end

nside = 8
res = Resolution(nside)
npix = nside2npix(nside)
l = zeros(npix)
b = zeros(npix)
for j in 1:npix
    theta, phi = pix2angRing(res, j)
    b[j] = 90 - rad2deg(theta)
    l[j] = rad2deg(phi)
end

size(out) == (npix, 3) || error("unexpected output size $(size(out))")
energies == [50.0, 100.0, 200.0] || error("unexpected energies $energies")
err_b = maximum(abs.(out[:, 1] .- b))
err_cos = maximum(abs.(out[:, 2] .- cos.(deg2rad.(l))))
err_sin = maximum(abs.(out[:, 3] .- sin.(deg2rad.(l))))
println("max |out - b| = $err_b, max |out - cos l| = $err_cos, max |out - sin l| = $err_sin")
err_b < 1e-9 || error("latitude plane not reproduced: $err_b")
err_cos < 2e-3 || error("cos(l) plane not reproduced: $err_cos")
err_sin < 2e-3 || error("sin(l) plane not reproduced across the seam: $err_sin")
println("template sampling ok")
