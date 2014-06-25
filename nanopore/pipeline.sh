#!/bin/bash

#Cleanup the old files
rm -rf ${2} ${3}
#Set the python path to just these local directories
export PYTHONPATH=:./:./submodules:${PYTHONPATH}
#Preferentially put the local binaries at the front of the path
export PATH=:./submodules/kmer:./submodules/sonLib/bin:./submodules/jobTree/bin:./submodules/bwa/:./submodules/lastz/src/:./submodules/last/src/:./submodules/last/scripts/:${PATH}
python2.7 ./nanopore/pipeline.py ${1} --jobTree ${2} --logInfo &> ${3}
