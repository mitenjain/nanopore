Code to analyse long read mappings, specifically tailored for Oxford Nanopore Reads.

### Requirements
* git
* python 2.7
* pysam

I'm pondering using virtualenv to eliminate the pysam dependency and allow other packages to be installed, but then you need to have virtualenv installed.

### Installation
To install the code run:

    git clone git://github.com/mitenjain/nanopore.git
    cd nanopore
    git pull
    git submodule update --init
    make

### Testing
To test the installation run:

    make test
    
This will run the demo sequences across the analyses in the test/ directory. The test sets mirrors the process of analysing your own data (see [Analysing your own data](https://github.com/mitenjain/nanopore/master/README.md#analysing-your-own-data)). 
    
### Updating the installation
To update a progressiveCactus installation, from the nanopore base directory type:

    git pull
    git submodule update --init
    make clean
    make

### Analysing your own data
The inputs are:
(1) One or more read files, in FASTA format. These should be placed in the directory *nanopore/readFastaFiles*.
(2) One or more reference genomes, in FASTA format. These should be placed in the directory *nanopore/referenceFastaFiles*.

For each possible pair of read file, reference genome and mapping algorithm an experiment directory will be created in the *nanopore/output* directory.

To run the pipeline from the nanopore base directory type:

    make run
    
To clean up an old run type:

    make clean

To see and control which mappers and analyses are being run edit lines 8 and 9 of the src/pipeline.py script.
    
