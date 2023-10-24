#!/bin/bash

JOB_ENVIRONMENT="fermi"
eval "$(conda shell.bash hook)"
conda activate $JOB_ENVIRONMENT

cd /lhome/ific/a/aamerio/github/MapGenData.jl
julia --threads $1 sub/create_artifacts_one_bin.jl