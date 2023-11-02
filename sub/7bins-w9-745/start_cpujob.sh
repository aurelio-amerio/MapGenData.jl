#!/bin/bash

JOB_ENVIRONMENT="fermi"
eval "$(conda shell.bash hook)"
conda activate $JOB_ENVIRONMENT

cd /lhome/ific/a/aamerio/github/MapGenData.jl
julia --threads $1 sub/7bins-w9-745/create_artifacts-w-9-745.jl