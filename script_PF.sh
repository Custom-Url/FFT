#!/bin/bash

# Trap Ctrl+C and kill all child processes
trap 'echo "Interrupted. Killing child processes..."; kill 0; exit 1' INT

# Set environment
source /mnt/data/dg765/FFT
export LD_LIBRARY_PATH=/mnt/data/dg765/FFT/amitex_fftp/libAmitex/lib:$LD_LIBRARY_PATH

# AMITEX executable and MPI
AMITEX="/mnt/data/dg765/FFT/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp"
MPIRUNalias='mpirun -np 18'

# Constant files
MATE="mate/mate_PF.xml"
LOAD="load/char.xml"
ALGO="algo/algo_Miehe2.xml"

# Loop through all generated VTK files in the nested directory structure
for vtk_file in micr/code_python/mesh/L*/h*/iUC_*.vtk; do
    # Extract spacing and volume fraction dirs
    spacing_dir=$(basename "$(dirname "$vtk_file")")        # h0.0002
    vf_dir=$(basename "$(dirname "$(dirname "$vtk_file")")") # L0.05_vf0.35

    resu_dir="resu/EXP/${vf_dir}/${spacing_dir}"
    mkdir -p "$resu_dir"

    echo "Running AMITEX for: $vtk_file"
    echo "Results in: $resu_dir"

    $MPIRUNalias $AMITEX -nm "$vtk_file" -m "$MATE" -c "$LOAD" -a "$ALGO" -s "$resu_dir"
done


