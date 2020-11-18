#!/bin/bash
NAME=${1:-$SLURM_JOB_NAME}
export TERACHEM_INP=$NAME.inp
export TERACHEM_EXE=tc.sh

mdrun -nt 1 -s $NAME.tpr