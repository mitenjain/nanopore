from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from jobTree.src.bioio import fastqRead, fastaRead, system, reverseComplement
from nanopore.analyses.utils import samIterator, getFastaDictionary, UniqueList
import pysam, os, itertools
from collections import Counter
from math import log

class IndelKmerAnalysis(AbstractAnalysis):
    """Runs kmer analysis"""

    def indelKmerFinder(self, aligned):
        r = UniqueList(); s = self.kmerSize+1
        for i in xrange(len(aligned)):
            r.add(aligned[i])
            if r[0] == None or (len(r) == s and r[self.kmerSize] == None) or (None not in r and len(r) == s):
                r.remove(last=False)
            elif None in r and len(r) == s:
                yield (r[0],r[self.kmerSize])
                r.remove(last=False)

    def countIndelKmers(self):
        sam = pysam.Samfile(self.samFile)
        refKmers, readKmers = Counter(), Counter()

        refDict = getFastaDictionary(self.referenceFastaFile)
        for x in refDict:
            refDict[x] = tuple(refDict[x])

        for record in samIterator(sam):
            refSeq = refDict[sam.getrname(record.rname)]
            readSeq = tuple(record.query)
            readAligned, refAligned = zip(*record.aligned_pairs)
            for start, end in self.indelKmerFinder(readAligned):
                if record.is_reverse:
                    readKmer = readSeq[start:end+1][::-1]
                else:
                    readKmer = readSeq[start:end+1]
                readKmers[readKmer] += 1
            for start, end in self.indelKmerFinder(refAligned):
                refKmers[refSeq[start:end+1]] += 1
        return (refKmers, readKmers)


    def analyzeCounts(self, refKmers, readKmers, name):
        refSize, readSize = sum(refKmers.values()), sum(readKmers.values())
        outf = open(os.path.join(self.outputDir, name + "kmer_counts.txt"), "w")
        outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tlogFoldChange\n")
        if refSize > 0 and readSize > 0:
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
        
            system("Rscript nanopore/analyses/kmer_analysis.R {} {} {}".format(os.path.join(self.outputDir, name + "kmer_counts.txt"), os.path.join(self.outputDir, name + "pval_kmer_counts.txt"), os.path.join(self.outputDir, name + "top_bot_sigkmer_counts.txt")))

    def run(self, kmerSize=5):
        AbstractAnalysis.run(self)
        self.kmerSize = kmerSize

        #analyze kmers around the boundaries of indels
        indelRefKmers, indelReadKmers = self.countIndelKmers()
        if len(indelRefKmers) > 0 and len(indelReadKmers) > 0:
            self.analyzeCounts(indelRefKmers, indelReadKmers, "indel_bases_")
        self.finish()