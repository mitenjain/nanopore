
all :
	cd submodules && make all

run : all
	./src/pipeline.sh ./ jobTree log.txt 

test : all
	./src/pipeline.sh tests testJobTree testLog.txt

clean :
	cd submodules && make clean
	rm -rf output jobTree log.txt tests/output testJobTree testLog.txt