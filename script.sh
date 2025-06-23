#!/bin/bash
#trap 'exit' INT


# ADJUST ENVIRONMENT VARIABLES if necessary ('module load xxx' on clusters)
# for example:
source /u/q/dg765/amitex_fftp/env_amitex.sh

# EXECUTABLE amitex_fftp
AMITEX="/u/q/dg765/amitex_fftp/libAmitex/src/user_Miehe2/amitex_fftp"

MPIRUNalias='mpirun -np 18'


SAMPLE="FM_vf0.6"
VXSIZ="h0.0001"
MICR="micr/"$SAMPLE"/"$VXSIZ".vtk"
MATE="mate/mate_2.xml"
LOAD="load/char.xml"
ALGO="algo/algo.xml"

mkdir "resu/"$SAMPLE
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s "resu/"$SAMPLE"/"$VXSIZ



