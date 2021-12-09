#!/bin/bash
# Load the necessary modules (software)

module purge
module load ncarenv/1.3
module load nvhpc/21.9
module load ncarcompilers/0.5.0
module load openmpi/4.1.1
module load netcdf/4.8.0
module load cuda/11.0.3
module list

# Export variables for use in the Makefile
export CUDA_ROOT_PATH="${NCAR_ROOT_CUDA}"
export NVHPC_ROOT_PATH="${NCAR_ROOT_NVHPC}/Linux_x86_64/20.11/compilers"

export PGI_ACC_TIME=1

./mkmf_preprocess
make clean
make
./preprocess

./mkmf_test_get_close_obs -mpi
make clean
make
