from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
from jobTree.src.bioio import system
from itertools import product

class UnmappedKmerAnalysis(AbstractUnmappedMetaAnalysis):
    """runs kmer analysis on all mapped/unmapped per read Type"""
    def run(self, kmer_size=5):
    	for readType in self.info.readTypes:
            unmapped = open(os.path.join(self.outputDir, readType + "_unmapped"), "w")
            mapped = open(os.path.join(self.outputDir, readType + "_mapped"), "w")
            for read in self.reads:
                if read.is_mapped and read.readType == readType:
                    mapped.write(">{}\n{}\n".format(read.name, read.seq))
                elif read.readType == readType:
                    unmapped.write(">{}\n{}\n".format(read.name, read.seq))
            unmapped.close(); mapped.close()
            system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_unmapped"), os.path.join(self.outputDir, readType + "_unmapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_mapped"), os.path.join(self.outputDir, readType + "_mapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_mapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_unmapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_" + str(kmer_size) + "kmer_Cmp.out")))
            system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(os.path.join(self.outputDir, readType + "_" + str(kmer_size)) + "kmer_Cmp.out"), os.path.join(self.outputDir, readType + "_top_kmers.tsv"), os.path.join(self.outputDir, readType + "_bot_kmers.tsv")))
