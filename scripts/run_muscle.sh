#!/bin/bash

while [[ $# > 1 ]]
do
key="$1"
shift

case $key in
    --template_sam)
    TEMPLATE_SAM="$1"
    shift
    ;;
    --twoD_sam)
    TWOD_SAM="$1"
    shift
    ;;
    --complement_sam)
    COMPLEMENT_SAM="$1"
    shift
    ;;
esac
done

if [ ! -e $COMPLEMENT_SAM ]; then
    echo "$COMPLEMENT_SAM does not exist"
    exit 1
fi
if [ ! -e $TEMPLATE_SAM ]; then
    echo "$TEMPLATE_SAM does not exist"
    exit 1
fi
if [ ! -e $TWOD_SAM ]; then
    echo "$TWOD_SAM does not exist"
    exit 1
fi

if [ -e ./jobTree ]; then
    rm -rf ./jobTree
fi

maxThreads=4
batchSystem=parasol
defaultJobMemory=8589934592

#Set the python path to just these local directories
export PYTHONPATH=:../:../submodules:${PYTHONPATH}
#Preferentially put the local binaries at the front of the path
export PATH=:../submodules/jobTree/bin:./submodules/muscle/:${PATH}:
python muscle_compare_2d/muscle_compare_2d.py $TEMPLATE_SAM $COMPLEMENT_SAM $TWOD_SAM --batchSystem=$batchSystem --jobTree=jobTree --defaultMemory=$defaultJobMemory --logInfo --logLevel INFO --stats &> log.txt

