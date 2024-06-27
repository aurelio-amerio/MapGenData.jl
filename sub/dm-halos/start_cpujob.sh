#!/bin/bash
cd /lhome/ific/a/aamerio/github/MapGenData.jl
julia --threads $1 sub/dm-halos/create_artifacts.jl