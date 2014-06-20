from nanopore.src.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system

class Bwa(AbstractMapper):
    def run(self):
        system("bwa index %s" % self.referenceFastaFile)
        system("bwa mem -x pacbio %s %s > %s" % (self.referenceFastaFile, self.readFastaFile, self.outputSamFile))