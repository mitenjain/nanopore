
all :
	echo Does nothing

run :
	rm -rf jobTree
	python src/pipeline.py ./ --logInfo > log.txt 2>&1

test :
	rm -rf testJobTree testLog.txt
	python src/pipeline.py tests --jobTree testJobTree --logInfo > testLog.txt 2>&1

clean :
	rm -rf output jobTree log.txt tests/output testJobTree testLog.txt