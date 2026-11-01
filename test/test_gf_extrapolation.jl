# Standalone check of the energy interpolation of the smoothed foreground template
# (interpolate_galactic_fg, used by the counts step).
# Run: julia +1.10 --project=<MapGenData.jl> test/test_gf_extrapolation.jl
# Each pixel follows a power law f = A E^-s through three planes; the interpolant must
# reproduce it inside the planes and, linearly in log-log, beyond the first and last
# plane. Interpolations.jl's (NoInterp(), Gridded(Linear())) with (Throw(), Line())
# holds the value flat outside the knots, which this check rejects.
using MapGenData
using StaticArrays
using Unitful

nside = 1
npix = 12 * nside^2
art = MapGenData.JLD2Artifact{1}("", nside, SVector(1.0), SVector(2.0))  # only nside is used
energies = [1.0e5, 3.0e5, 8.0e5]                                          # MeV
A = [1.0e-6 * (1 + p) for p in 1:npix]
s = [2.0 + 0.1 * p for p in 1:npix]
model = [A[p] * E^-s[p] for p in 1:npix, E in energies]                   # (npix, nE)

f, E_out = MapGenData.interpolate_galactic_fg(model, energies, art)
ustrip.(u"MeV", E_out) == energies || error("unexpected energies $E_out")

maxrel = 0.0
for p in 1:npix, E in (5.0e4, 1.0e5, 2.0e5, 8.0e5, 9.0e5, 1.0e6)       # below, inside, at and above the planes
    got = ustrip(u"MeV^-1", f(p, E * u"MeV"))
    want = A[p] * E^-s[p]
    global maxrel = max(maxrel, abs(got - want) / want)
end
println("max relative error of the log-log interpolation (inside and outside the planes) = $maxrel")
maxrel < 1e-12 || error("template not interpolated/extrapolated linearly in log-log: max rel error $maxrel")
println("gf extrapolation ok")
