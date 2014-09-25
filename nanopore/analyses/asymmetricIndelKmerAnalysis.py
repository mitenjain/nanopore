from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from jobTree.src.bioio import fastqRead, fastaRead, system, reverseComplement
from nanopore.analyses.utils import samIterator, getFastaDictionary
import pysam, os, itertools
from collections import Counter
from math import log


class AsymmetricIndelKmerAnalysis(AbstractAnalysis):
    """Runs kmer analysis"""
    def countIndelKmers(self):
        sam = pysam.Samfile(self.samFile)
        refKmers, readKmers = Counter(), Counter()

        refDict = getFastaDictionary(self.referenceFastaFile)
        for x in refDict:
            refDict[x] = tuple(refDict[x])

        aRs = list()
        for record in samIterator(sam):
            aRs.append([record.aligned_pairs, tuple(record.query), sam.getrname(record.rname), record.is_reverse])
            #trying to be fast while not running out of RAM
            if len(aRs) % 5000 == 0:
                for aP, readSeq, name, is_reverse in aRs:
                    refSeq = refDict[name]
                    refKmer, readKmer, cutoff = list(), list(), 0
                    for i in xrange(len(aP)):
                        readPos, refPos = aP[i]
                        refKmer.append(refPos); readKmer.append(readPos)
                        #if we hit kmerSize bases without indels, start saving the cutoff
                        #this is faster than making a new list every time
                        if len(refKmer) == 5 and None not in refKmer and None not in readKmer:
                            cutoff += 1
                        #we are exiting a window of indels
                        elif None not in refKmer[-5:] and None not in readKmer[-5:] and len(refKmer) > 5:
                            refKmer = [x for x in refKmer[cutoff:] if not x == None]
                            readKmer = [x for x in readKmer[cutoff:] if not x == None]
                            for i in refKmer[5:]:
                                refKmers[refSeq[i-5:i]] += 1
                            for i in readKmer[5:]:
                                seq = readSeq[i-5:i]
                                if is_reverse:
                                    readKmers[seq[::-1]] += 1
                                else:
                                    readKmers[seq] += 1
                            refKmer, readKmer, cutoff = list(), list(), 0
            aRs = list()
        #need to get the last 5000 records
        for aP, readSeq, name, is_reverse in aRs:
            refSeq = refDict[name]
            refKmer, readKmer, cutoff = list(), list(), 0
            readPositions, refPositions = zip(*aP)
            for i in xrange(self.kmerSize, len(aP)):
                if None in readPositions[i-self.kmerSize:i] and None not in refPositions[i-self.kmerSize:i]:
                    refKmers[refSeq[i-self.kmerSize:i]] += 1
                if None not in readPositions[i-self.kmerSize:i] and None in refPositions[i-self.kmerSize:i]:
                    seq = readSeq[i-self.kmerSize:i]
                    if is_reverse:
                        readKmers[seq[::-1]] += 1
                    else:
                        readKmers[seq] += 1
        return (refKmers, readKmers) 

    def analyzeCounts(self, refKmers, readKmers, name):
        refSize, readSize = sum(refKmers.values()), sum(readKmers.values())
        outf = open(os.path.join(self.getLocalTempDir(), name + "kmer_counts.txt"), "w")
        outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tlogFoldChange\n")
        
        for kmer in itertools.product("ATGC", repeat=5):
            refFraction, readFraction = 1.0 * refKmers[kmer] / refSize, 1.0 * readKmers[kmer] / readSize
            if refFraction == 0:
                foldChange = "-Inf"
            elif readFraction == 0:
                foldChange = "Inf"
            else:
                foldChange = -log(readFraction / refFraction)
            outf.write("\t".join(map(str,["".join(kmer), refKmers[kmer], refFraction, readKmers[kmer], readFraction, foldChange]))+"\n")
        outf.close()
        
        system("Rscript nanopore/analyses/kmer_analysis.R {} {} {}".format(os.path.join(self.getLocalTempDir(), name + "kmer_counts.txt"), os.path.join(self.outputDir, name + "kmer_counts.txt"), os.path.join(self.outputDir, name + "top_bot_sigkmer_counts.txt")))

    def run(self, kmerSize=5):
        AbstractAnalysis.run(self)
        self.kmerSize = kmerSize

        #analyze kmers around the boundaries of indels
        indelRefKmers, indelReadKmers = self.countIndelKmers()
        if len(indelRefKmers) > 0 and len(indelReadKmers) > 0:
            self.analyzeCounts(indelRefKmers, indelReadKmers, "asymmetric_indel_bases_")
        self.finish()