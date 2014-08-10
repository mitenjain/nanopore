maxThreads = 4
batchSystem = singleMachine
defaultJobMemory=8589934592

all :
	cd submodules && make all

run : all
	./nanopore/pipeline.sh ./ jobTree log.txt ${maxThreads} ${batchSystem} ${defaultJobMemory}

test : all
	./nanopore/pipeline.sh tests testJobTree testLog.txt ${maxThreads} ${batchSystem} ${defaultJobMemory}

clean :
	cd submodules && make clean
	rm -rf output jobTree log.txt tests/output testJobTree testLog.txt
