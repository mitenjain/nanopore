maxThreads = 4
batchSystem = singleMachine

all :
	cd submodules && make all

run : all
	./nanopore/pipeline.sh ./ jobTree log.txt ${maxThreads} ${batchSystem}

test : all
	./nanopore/pipeline.sh tests testJobTree testLog.txt ${maxThreads} ${batchSystem}

clean :
	cd submodules && make clean
	rm -rf output jobTree log.txt tests/output testJobTree testLog.txt