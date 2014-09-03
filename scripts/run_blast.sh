#!/bin/bash

maxThreads=20
batchSystem=parasol
defaultJobMemory=8589934592

if [ -e ./jobTree ]; then
    rm -rf ./jobTree
fi


export BLASTDB=:/hive/users/ifiddes/blastdb
#Set the python path to just these local directories
export PYTHONPATH=:../:../submodules:${PYTHONPATH}
#Preferentially put the local binaries at the front of the path
export PATH=:../submodules/jobTree/bin:./submodules/muscle/:${PATH}:
python blast_combined/blast_combined.py --batchSystem=$batchSystem --jobTree=jobTree --defaultMemory=$defaultJobMemory --logInfo --logLevel INFO --stats &> log.txt
