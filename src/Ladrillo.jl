"""
    Ladrillo

A calibrated sea-level emulator: MimiBRICK v2.0.0 with a three-reservoir glacier component and a two-basin Greenland
component, calibrated jointly to observed component sea-level budgets (posterior = calibration L27, shared by v1.0 and v1.1).

Quick start:
```julia
using Ladrillo
gmst, ohc = Ladrillo.mean_forcing("ssp245")
lm = LadrilloModel(gmst, ohc)
post = Ladrillo.posterior()
run_draw!(lm, post[1, :])
sea_level(lm, :total)          # cm rel. 1995-2014, 1850-2300
```
"""
module Ladrillo

using CSV, DataFrames, Mimi, MimiBRICK, Random, Statistics

export LadrilloModel, apply_draw!, run_draw!, sea_level, posterior, mean_forcing, forcing_cube, splice

const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
datapath(parts...) = joinpath(DATA_DIR, parts...)

"""Lazily filled values read from data/ (kept out of precompilation so a data update cannot be shadowed by a stale
compile cache)."""
const _CACHE = Dict{Symbol,Any}()

const GLACIER_BLOCKS = ("R19", "SLOWG", "FASTG")

include("components/glaciers.jl")
include("components/greenland.jl")
include("landwater.jl")
include("posterior.jl")
include("model.jl")
include("forcing.jl")
include("targets.jl")
include("brick20.jl")

end # module
