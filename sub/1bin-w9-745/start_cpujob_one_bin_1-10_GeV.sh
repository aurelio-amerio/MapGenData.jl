#!/bin/bash

JOB_ENVIRONMENT="fermi"
eval "$(conda shell.bash hook)"
conda activate $JOB_ENVIRONMENT

cd /lhome/ific/a/aamerio/github/MapGenData.jl
julia --threads $1 sub/1bin-w9-745/create_artifacts_one_bin_1-10_GeV.jl