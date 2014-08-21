from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class Bwa(AbstractMapper):
    def run(self, args=""):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because BWA builds these crufty index files, copy to a temporary directory
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile))
        system("bwa index %s" % localReferenceFastaFile)
        system("bwa mem %s %s %s > %s" % (args, localReferenceFastaFile, self.readFastqFile, self.outputSamFile))

class BwaChain(Bwa):
    def run(self):
        Bwa.run(self)
        self.chainSamFile()

class BwaRealign(Bwa):
    def run(self):
        Bwa.run(self)
        self.realignSamFile()

class BwaRealignEm(Bwa):
    def run(self):
        Bwa.run(self)
        self.realignSamFile()