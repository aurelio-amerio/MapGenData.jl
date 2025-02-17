#!/bin/bash

# JOB_ENVIRONMENT="fermi"
# eval "$(conda shell.bash hook)"
# conda activate $JOB_ENVIRONMENT

cpus=$1
threads=$((2*cpus))

cd /lhome/ific/a/aamerio/github/MapGenData.jl
julia --threads $threads sub/7bins/create_artifacts.jl