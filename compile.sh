#!/bin/bash

# Exit immediately if a command fails
set -e

# 1. Load the correct modules
module load mpi
module load gnu8/8.3.0

# 2. Source clik environment
source $HOME/scope/mcmc/code/plc_3.0/plc-3.1/bin/clik_profile.sh

# 3. Define paths for clarity
CAMB_LIB_DIR="$HOME/scope/mcmc/CAMB/fortran/Releaselib"
LIBS_DIR="$HOME/scope/libs/lib"
CLIK_LIB_DIR="$HOME/scope/mcmc/code/plc_3.0/plc-3.1/lib"

echo "Compiling mcmc1..."

# 4. Compile Command
# Note: We pass the full path to camblib.so because it is not named 'libcamblib.so',
# so -lcamblib would not find it.
mpicc -std=gnu99 -O3 -fopenmp -o mcmc1 mcmc.c \
      -I$HOME/scope/libs/include \
      -I$HOME/scope/mcmc/code/plc_3.0/plc-3.1/src/minipmc \
      $CAMB_LIB_DIR/camblib.so \
      -L$LIBS_DIR \
      -L$CLIK_LIB_DIR \
      -L/opt/ohpc/pub/compiler/gcc/8.3.0/lib64 \
      -lclik -lgsl -lgslcblas -lcfitsio -llapack -lblas \
      -lgfortran \
      -lpthread -lm -lrt \
      -Wl,-rpath,$CAMB_LIB_DIR \
      -Wl,-rpath,$LIBS_DIR \
      -Wl,-rpath,$CLIK_LIB_DIR

echo "--- Compile successful! ---"