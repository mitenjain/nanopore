#!/bin/bash

PYTHON=`which python`
rm -rf ${2} ${3}
export PYTHONPATH=:../:./subModules
export PATH=:../nanopore/subModules/sonLib/bin:../nanopore/subModules/jobTree/bin:${PATH}
$PYTHON src/pipeline.py ${1} --jobTree ${2} --logInfo &> ${3}