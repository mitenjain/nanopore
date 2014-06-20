# order is important, libraries first
modules = sonLib jobTree bwa lastz last
.PHONY: all %.all clean %.clean

all : ${modules:%=all.%} 

all.%:
	cd $* && make

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean