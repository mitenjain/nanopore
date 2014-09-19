from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
from jobTree.src.bioio import system
import itertools
from collections import Counter
from math import log

class UnmappedKmerAnalysis(AbstractUnmappedMetaAnalysis):
    """Calculates kmer statistics for all reads (in all samples) not mapped by any mapper
    """

    def countKmers(self, seq):
        kmers = Counter()
        for i in xrange(self.kmerSize, len(seq)):
            if "N" not in seq[ i - self.kmerSize : i ]:
                kmers[ seq[ i - self.kmerSize : i ] ] += 1
        return kmers

    def run(self, kmerSize=5):
        self.kmerSize = kmerSize
        for readType in self.readTypes:
            mappedKmers, unmappedKmers = Counter(), Counter()
            for read in self.reads:
                if read.readType == readType and read.is_mapped:
                    mappedKmers += self.countKmers(read.seq)
                elif read.readType == readType:
                    unmappedKmers += self.countKmers(read.seq)

            mappedSize, unmappedSize = sum(mappedKmers.values()), sum(unmappedKmers.values())
            outf = open(os.path.join(self.getLocalTempDir(), readType + "_kmer_counts.txt"), "w")
            outf.write("kmer\tmappableCount\tmappableFraction\tunmappableCount\tunmappableFraction\tlogFoldChange\n")
            for kmer in itertools.product("ATGC",repeat=5):
                kmer = "".join(kmer)
                mappedFraction, unmappedFraction = 1.0 * mappedKmers[kmer] / mappedSize, 1.0 * unmappedKmers[kmer] / unmappedSize
                foldChange = -log(mappedFraction / unmappedFraction)
                outf.write("\t".join(map(str,[kmer, mappedKmers[kmer], mappedFraction, unmappedKmers[kmer], unmappedFraction, foldChange]))+"\n")
            outf.close()

            system("Rscript nanopore/analyses/kmer_analysis.R {} {} {}".format(os.path.join(self.getLocalTempDir(), readType + "_kmer_counts.txt"), os.path.join(self.outputDir, readType + "_unmapped_kmer_counts.txt"), os.path.join(self.outputDir, readType + "_unmapped_top_bot_sigkmer_counts.txt")))
