#!/bin/bash

#Cleanup the old files
rm -rf ${2} ${3}
#Set the python path to just these local directories
export PYTHONPATH=:./:./submodules:${PYTHONPATH}
#Preferentially put the local binaries at the front of the path
export PATH=:./submodules/sonLib/bin:./submodules/cactus/bin:./submodules/jobTree/bin:./submodules/bwa/:./submodules/lastz/src/:./submodules/last/src/:./submodules/last/scripts/:./submodules/fastqc/:./submodules/qualimap:./submodules/blasrbinary/:./submodules/samtools-0.1.19:./submodules/samtools-0.1.19/bcftools:./submodules/consensus/:${PATH}:
python ./nanopore/pipeline.py ${1} --jobTree ${2} --logInfo --maxThreads=${4} --batchSystem=${5} --defaultMemory=${6} --logLevel INFO --stats &> ${3}
