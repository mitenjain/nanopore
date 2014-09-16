from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from jobTree.src.bioio import fastqRead, fastaRead, system, reverseComplement
from nanopore.analyses.utils import AlignedPair
import pysam, os, itertools
from collections import Counter
from math import log


class KmerAnalysis(AbstractAnalysis):
"""Runs kmer analysis"""
    def countIndelKmers(self):
        sam = pysam.Samfile(self.sam)
        refKmers, readKmers = Counter(), Counter()
        for aR in samIterator(sam):
            #Ref name
            refSeqName = sam.getrname(aR.rname)
            #Sequences
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = aR.query
            #temp lists of kmers
            refKmer, readKmer = list(), list()
            #flag that says if we are currently in an indel
            indelFlag = False
            #iterate over aligned pairs of positions
            for aP in AlignedPair.iterator(aR, refSeq, readSeq):
                #add the base at each position to the list
                refKmer.append(aP.getRefBase()); readKmer.append(aP.getReadBase())
                #if this position has a indel, we want to save this kmer as soon as it hits kmerSize
                if aP.getPrecedingReadInsertionLength != 0 or aP.getPrecedingReadDeletionLength != 0:
                    indelFlag = True
                if indelFlag is True and len(refKmer) == self.kmerSize:
                    if aP.isReversed:
                        refKmers["".join(refKmer)] += 1; readKmers[reverseComplement("".join(readKmer))] += 1
                    else:
                        refKmers["".join(refKmer)] += 1; readKmers["".join(readKmer)] += 1
                    indelFlag = False
        return (refKmers, readKmers)

    def countKmers(self):
        refKmers, readKmers = Counter(), Counter()

        for name, seq in fastaRead(self.referenceFastaFile):
            for i in xrange(kmerSize, len(seq)):
                if "N" not in seq[ i - kmerSize : i ]:
                    refKmers[ seq[ i - kmerSize : i ] ] += 1

        for name, seq, qual in fastqRead(self.readFastqFile):
            for i in xrange(kmerSize, len(seq)):
                if "N" not in seq[ i - kmerSize : i ]:
                    readKmers[ seq[ i - kmerSize : i ] ] += 1
        return (refKmers, readKmers)

    def analyzeCounts(self, refKmers, readKmers, name):
        refSize, readSize = sum(refKmers.values()), sum(readKmers.values())
        outf = open(os.path.join(self.getLocalTempDir(), name + "kmer_counts.txt"), "w")
        outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tfoldChange\n")
        for kmer in itertools.product("ATGC",repeat=5):
            kmer = "".join(kmer)
            refFraction, readFraction = 1.0 * refKmers[kmer] / refSize, 1.0 * readKmers[kmer] / readSize
            foldChange = -log(readFraction / refFraction)
            outf.write("\t".join(map(str,[kmer, refKmers[kmer], refFraction, readKmers[kmer], readFraction, foldChange]))+"\n")
        outf.close()

        system("Rscript nanopore/analyses/kmer_analysis.R {} {}".format(os.path.join(self.getLocalTempDir(), name + "kmer_counts.txt"), os.path.join(self.outputDir, name + "kmer_counts.txt")))

    def run(self, kmerSize=5):
        AbstractAnalysis.run(self)
        self.kmerSize = kmerSize
        
        #analyze kmers across both files
        refKmers, readKmers = countKmers()
        analyzeCounts(refKmers, readKmers, "all_bases_")

        #analyze kmers around the boundaries of indels
        refKmers, readKmers = countIndelKmers()
        analyzeCounts(refKmers, readKmers, "indel_bases_")