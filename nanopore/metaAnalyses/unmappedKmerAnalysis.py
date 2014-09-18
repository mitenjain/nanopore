from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
from jobTree.src.bioio import system
from itertools import product
from collections import Counter

class UnmappedKmerAnalysis(AbstractUnmappedMetaAnalysis):
    """Calculates kmer statistics for all reads (in all samples) not mapped by any mapper
    """

    def countKmers(self, seq):
        kmers = Counter()
        for i in xrange(kmerSize, len(seq)):
            if "N" not in seq[ i - kmerSize : i ]:
                kmers[ seq[ i - kmerSize : i ] ] += 1
        return kmers

    def run(self, kmer_size=5):
        for readType in self.readTypes:
            mappedKmers, unmappedKmers = Counter(), Counter()
            for read in self.reads:
                if read.readType == readType and read.is_mapped:
                    mappedKmers += self.countKmers(read.seq)
                elif read.readType == readType:
                    unmappedKmers += self.countKmers(read.seq)

            mappedSize, unmappedSize = sum(mappedKmers.values()), sum(unmappedKmers.values())
            outf = open(os.path.join(self.getLocalTempDir(), "kmer_counts.txt"), "w")
            outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tfoldChange\n")
            for kmer in itertools.product("ATGC",repeat=5):
                kmer = "".join(kmer)
                refFraction, readFraction = 1.0 * refKmers[kmer] / refSize, 1.0 * readKmers[kmer] / readSize
                foldChange = -log(readFraction / refFraction)
                outf.write("\t".join(map(str,[kmer, refKmers[kmer], refFraction, readKmers[kmer], readFraction, foldChange]))+"\n")
            outf.close()

        system("Rscript nanopore/analyses/kmer_analysis.R {} {} {}".format(os.path.join(self.getLocalTempDir(), "kmer_counts.txt"), os.path.join(self.outputDir, readType + "kmer_counts.txt"), os.path.join(self.outputDir, name + "top_bot_sigkmer_counts.txt")))
