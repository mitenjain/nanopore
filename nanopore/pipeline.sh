#!/bin/bash

#Cleanup the old files
rm -rf ${2} ${3}
#Set the python path to just these local directories
export PYTHONPATH=:../:./submodules
#Preferentially put the local binaries at the front of the path
export PATH=:./submodules/sonLib/bin:./submodules/jobTree/bin:./submodules/bwa/:./submodules/lastz/src/:${PATH}
python ./nanopore/pipeline.py ${1} --jobTree ${2} --logInfo &> ${3}