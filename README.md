Code to analyse long read mappings, specifically tailored for Oxford Nanopore Reads.

### Requirements
* git
* python 2.7
* pysam
* R 2.15.1 or newer
* Lattice package for R (http://cran.r-project.org/web/packages/lattice/index.html)

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
(1) One or more read files, in FASTQ format. These should be placed in the directory *nanopore/readFastaFiles/<X>/*.
Where <X> is the readType being analyzed. Generally, this is in the format 2D/, complement/, and template/.
(2) One or more reference genomes, in FASTA format. These should be placed in the directory *nanopore/referenceFastaFiles*.

For each possible pair of read file, reference genome and mapping algorithm an experiment directory will be created in the *nanopore/output* directory.

To run the pipeline from the nanopore base directory type:

    make run
    
To clean up an old run type:

    make clean

To see and control which mappers and analyses are being run edit lines 37 and 66 of the src/pipeline.py script.

###Setting up BLAST
Additionally, if you would like to BLAST the unmapped reads against NT (NCBI nucleotide database) you need to have BLAST installed with the BLASTDB environmental variable set to wherever you have stored the NT database. In addition, you want the taxonomy database so that the hits make sense.

You can download blastn from here:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
You will want to move blastn to somewhere that is in your path.
Second, set the environmental variable BLASTDB:
```export BLASTDB=:/path/to/blast/db/``` (you will want this in your bashrc)
Finally, download and untar all of the databases to your BLASTDB path:
ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.**.tar.gz
ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz

For all of the (19) nt files, you have to untar separately. Thanks NCBI. Dirty solution:
```for i in `echo $BLASTDB | cut -d ":" -f 2`*.tar.gz; do tar zxvf $i; done``` 

###Scripts
There is also a scripts directory where scripts can be dropped to analyze pipeline results outside of the pipeline. The jobTree can still be used, if you design your controlling shell script to change the paths as shown in the currently only external script, run_muscle.sh. These scripts will often depend on the directory structure of the pipeline to find sam files, xml files, etc.

###Current scripts
`run_muscle.sh` - has three inputs: `--template_sam`, `--twoD_sam`, `--complement_sam`. Takes these three sam files and looks for all read and reference files that these reads came from, and determines which reads were mappable as 2D but not mappable as template/complement. Then, the region of the reference where the 2D read aligned is extracted and MUSCLE is ran to try and get a global alignment between the template and complement and the corresponding 2D aligned region.
