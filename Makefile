
all : 
	rm -rf jobTree && python src/pipeline.py ./ &2> log.txt

test : 
	rm -rf testJobTree && python src/pipeline.py ./test --jobTree ./testJobTree pwd &2> testLog.txt
	
clean :
	rm -rf output jobTree log.txt test/output testJobTree testLog.txt