from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class Bwa(AbstractMapper):
    def run(self):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because BWA builds these crufty index files, copy to a temporary directory
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile))
        system("bwa index %s" % localReferenceFastaFile)
        system("bwa mem -x pacbio %s %s > %s" % (localReferenceFastaFile, self.readFastaFile, self.outputSamFile))