#!/bin/bash
#
# This script verifies submission structurefrom the COMS30006 Lattice Boltzmann
# coursework. This script should be run from the directory containing the files
# to be submitted.
#
# This script will unload all modules and source 'env.sh' (if present). It will
# then run `make` and check that there is an executable with the correct name.

set -e

EXE=d2q9-bgk
ENV=env.sh

module list |& tail -n +2

echo "Unloading all modules"
module purge &>/dev/null

echo "Loading default cuda module"
module load libs/cuda/10.0-gcc-5.4.0-2.26

if [ -r "$ENV" ]; then
    echo "Sourcing $ENV"
    source "$ENV"
else
    echo "No $ENV present, skipping"
fi

module list

echo "Cleaning old executable"
rm -f "$EXE"

if [ ! -r Makefile ]; then
    echo -e "\nERROR: Makefile not found."
    exit 1
fi

echo 'Running `make`:'
if ! make -B; then
    echo -e "\nERROR: Build failed. Are you missing any modules from $ENV?"
    exit 11
elif [ ! -r "$EXE" ]; then
    echo -e "\nERROR: Executable '$EXE' is not present."
    echo    "If your executable name is different, please change it to '$EXE' in your Makefile."
    exit 12
else
    make -s clean
    cat <<-EOM

	Submission check passed. 

	Please ensure that you submit all source files in the submission directory:
	    - Makefile
	    - d2q9-bgk.c
	    - env.sh
	    - Any other files needed to build that you've added, e.g. OpenCL kernel files.

	Please also submit your report with the filename 'report.pdf'.

	Note: Your code has NOT been run and this does NOT mean that the results validate. You should check correctness separately.
	EOM
fi

