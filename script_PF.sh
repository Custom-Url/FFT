#!/bin/bash
#trap 'exit' INT


# ADJUST ENVIRONMENT VARIABLES if necessary ('module load xxx' on clusters)
# for example:
source /mnt/data/dg765/FFT

export LD_LIBRARY_PATH=/mnt/data/dg765/FFT/amitex_fftp/libAmitex/lib:$LD_LIBRARY_PATH


# EXECUTABLE amitex_fftp
AMITEX="/mnt/data/dg765/FFT/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp"

MPIRUNalias='mpirun -np 18'


SAMPLE="FM_vf0.6"
VXSIZ="h0.0001"
MICR="micr/"$SAMPLE"/"$VXSIZ".vtk"
MATE="mate/mate_PF.xml"
LOAD="load/char.xml"
ALGO="algo/algo_Miehe2.xml"

mkdir "resu/"$SAMPLE
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s "resu/"$SAMPLE"/"$VXSIZ



