from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
import pysam
from itertools import izip

class UnmappedKmer(AbstractMetaAnalysis):
    """Calculates kmer statistics for all reads (in all samples) not mapped by any mapper"""
    def run(self, kmer_size=5):
        unmapped_by_ref = dict()
        for readFastqFile, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            #make sets of sequences from both mapped and all reads
            mappedSequences = {x.seq for x in samIterator(pysam.Samfile(os.path.join(resultsDir, "mapping.sam")))}
            readSequences = {x[1] for x in fastqRead(open(readFastqFile))}
            #add readSequences - mappedSequences to unmapped
            if referenceFastaFile not in unmapped_by_ref:
                unmapped_by_ref[referenceFastaFile] = set()
            unmapped_by_ref[referenceFastaFile] = unmapped_by_ref[referenceFastaFile].union(readSequences.difference(mappedSequences))

        #first need to make reads into fasta format because I am not recoding Karen's perl
        for referenceFastaFile, unmapped in unmapped_by_ref.iteritems():
            outf = open(os.path.join(self.getLocalTempDir(), "tmp.fasta"), "w")
            for seq, i in izip(unmapped, xrange(len(unmapped))):
                outf.write(">{}\n{}\n".format(i, seq))
            outf.close()
            #now we run kmer analysis
            ref_base = os.path.basename(referenceFastaFile).split(".")[0]
            system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "tmp.fasta"), os.path.join(self.getLocalTempDir(), ref_base + "_read_" + str(kmer_size) + "mer"), str(kmer_size)))
            system("nanopore/analyses/kmer.pl {} {} {}".format(referenceFastaFile, os.path.join(self.getLocalTempDir(), ref_base + "_ref_" + str(kmer_size) + "mer"), str(kmer_size)))
            system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), ref_base + "_ref_" + str(kmer_size) + "mer"), os.path.join(self.getLocalTempDir(), ref_base + "_read_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, ref_base + "_" + str(kmer_size) + "kmer_Cmp.out")))
            system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(self.outputDir, ref_base + "_" + str(kmer_size) + "kmer_Cmp.out"), os.path.join(self.outputDir, "top_kmers.tsv"), os.path.join(self.outputDir, "bot_kmers.tsv")))

