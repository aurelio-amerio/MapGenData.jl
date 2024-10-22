using ArtifactUtils
using Base.Filesystem: cp
using Base.Threads 
using CodecZlib
using Downloads
using FITSIO
using HDF5
using Healpix
using Interpolations
using JLD2
using LegendrePolynomials
using Memoize
using MultiQuad
using Pkg
using ProgressMeter
using QuadGK
using Scratch
using StaticArrays
using Unitful
using Unitful: Energy, EnergyFreeUnits

import Base: convert, Matrix