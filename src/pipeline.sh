#!/bin/bash

PYTHON=`which python`
rm -rf ${2} ${3}
export PYTHONPATH=:../:./subModules
export PATH=:./subModules/sonLib/bin:./subModules/jobTree/bin:./subModules/bwa/:./subModules/lastz/src/:${PATH}
$PYTHON src/pipeline.py ${1} --jobTree ${2} --logInfo &> ${3}