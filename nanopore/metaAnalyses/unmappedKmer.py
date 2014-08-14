from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
import pysam

class UnmappedKmer(AbstractMetaAnalysis):
    """Calculates kmer statistics for all reads (in all samples) not mapped by any mapper
    This class works almost identically to UnmappedBlastKmer but without the blasting part
    And is supposed to be run on machines where blast is not installed (or if you don't want
    the pipeline to takes days)"""
    def run(self, kmer_size=5):
        #build set of (read, fastq) of all mapped reads
        mappedReads, unmappedReads = dict(), dict()
        for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            if readType not in mappedReads:
                mappedReads[readType] = set()
            for record in samIterator(pysam.Samfile(os.path.join(resultsDir, "mapping.sam"))):
                if not record.is_unmapped: #some aligners save unmapped reads in their samfiles (bwa)
                    mappedReads[readType].add((record.qname, readFastqFile))
        
        #loop again and save only reads that did not map with any mapper
        for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            if readType not in unmappedReads:
                unmappedReads[readType] = set()
            for name, seq, qual in fastqRead(readFastqFile):
                if (name, readFastqFile) not in mappedReads:
                    unmappedReads[readType].add((name, seq))

        #build fastas of mapped and unmapped reads for kmer analysis
        for readType in unmappedReads:
            if len(unmappedReads[readType]) >= 1:
                #first we build a fasta of unmapped reads
                outf = open(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), "w")
                for name, readFastqFile, seq in unmappedReads[readType]:
                    outf.write(">{}\n{}\n".format(name, seq))
                outf.close()
                #then we build a fasta of mapped reads, iterating over the fastq files so we get the full read:
                outf = open(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), "w")
                for readFastqFile, thisReadType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
                    if thisReadType == readType:
                        for name, seq, qual in fastqRead(readFastqFile):
                            if (name, readFastqFile) in mappedReads[readType]:
                                outf.write(">{}\n{}\n".format(name, seq))
                outf.close()
                system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), str(kmer_size)))
                system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), str(kmer_size)))
                system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_" + str(kmer_size) + "kmer_Cmp.out")))
                system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(os.path.join(self.outputDir, readType + "_" + str(kmer_size)) + "kmer_Cmp.out"), os.path.join(self.outputDir, readType + "_top_kmers.tsv"), os.path.join(self.outputDir, readType + "_bot_kmers.tsv")))
          