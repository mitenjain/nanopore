from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from jobTree.src.bioio import fastqRead, fastaRead, system, reverseComplement
from nanopore.analyses.utils import samIterator
import pysam, os, itertools
from collections import Counter
from math import log


class AsymmetricIndelKmerAnalysis(AbstractAnalysis):
    """Runs kmer analysis"""
    def countIndelKmers(self):
        sam = pysam.Samfile(self.samFile)
        refKmers, readKmers = Counter(), Counter()
        
        for aR in samIterator(sam):
            refSeq = None
            for name, seq in fastaRead(self.referenceFastaFile):
                if name == sam.getrname(aR.rname):
                    refSeq = seq
            if refSeq is None:
                continue

            readSeq = aR.query
            #temp lists of kmers
            refKmer, readKmer = list(), list()
            for i in xrange(len(aR.aligned_pairs)):
                readPos, refPos = aR.aligned_pairs[i]
                refKmer.append(refPos); readKmer.append(readPos)
                
                #if we hit kmerSize bases without indels, start tossing out positions
                if len(refKmer) == self.kmerSize and None not in refKmer and None not in readKmer:
                    refKmer = refKmer[1:]; readKmer = readKmer[1:]
                
                #we are exiting a window of indels
                elif None not in refKmer[-self.kmerSize:] and None not in readKmer[-self.kmerSize:] and len(refKmer) > self.kmerSize:
                    refKmer = [x for x in refKmer if not x == None]
                    readKmer = [x for x in readKmer if not x == None]
                    for i in refKmer[self.kmerSize:]:
                        refKmers[refSeq[i-self.kmerSize:i]] += 1
                    for i in readKmer[self.kmerSize:]:
                        seq = readSeq[i-self.kmerSize:i]
                        if aR.is_reverse:
                            readKmers[reverseComplement(seq)] += 1
                        else:
                            readKmers[seq] += 1
                    refKmer, readKmer = list(), list()
        return (refKmers, readKmers)    

    def analyzeCounts(self, refKmers, readKmers, name):
        refSize, readSize = sum(refKmers.values()), sum(readKmers.values())
        outf = open(os.path.join(self.getLocalTempDir(), name + "kmer_counts.txt"), "w")
        outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tlogFoldChange\n")
        
        for kmer in itertools.product("ATGC", repeat=5):
            kmer = "".join(kmer)
            refFraction, readFraction = 1.0 * refKmers[kmer] / refSize, 1.0 * readKmers[kmer] / readSize
            if refFraction == 0:
                foldChange = "-Inf"
            elif readFraction == 0:
                foldChange = "Inf"
            else:
                foldChange = -log(readFraction / refFraction)
            outf.write("\t".join(map(str,[kmer, refKmers[kmer], refFraction, readKmers[kmer], readFraction, foldChange]))+"\n")
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