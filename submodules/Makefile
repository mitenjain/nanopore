# order is important, libraries first
modules = sonLib jobTree bwa lastz last
.PHONY: all %.all clean %.clean

export HDF5INCLUDEDIR=$(PWD)/hdf5Build/include
export HDF5LIBDIR=$(PWD)/hdf5Build/lib

all : ${modules:%=all.%} 

all.%:
	cd $* && make

#This is needed for blasr, but currently blasr won't compile on either hgwdev or os x
hdf5Rule :
	cd hdf5 && ./configure --prefix=$(PWD)/hdf5Build --enable-cxx && make && make install

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM
	cd hdf5 && make clean
	rm -rf hdf5Build

clean.%:
	cd $* && make clean