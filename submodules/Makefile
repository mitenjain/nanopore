# order is important, libraries first
modules = sonLib jobTree bwa lastz last cPecan samtools-0.1.19 
#hdf5 blasr
.PHONY: all %.all clean %.clean

export HDF5INCLUDEDIR=$(PWD)/hdf5/hdf5/include
export HDF5LIBDIR=$(PWD)/hdf5/hdf5/lib

all : ${modules:%=all.%} 

all.%:
	cd $* && make

#This is needed for blasr, but currently blasr won't compile on either hgwdev or os x
hdf5Rule :
	cd hdf5 && ./configure --enable-cxx && make && make install

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM
	cd hdf5 && make clean

clean.%:
	cd $* && make clean
