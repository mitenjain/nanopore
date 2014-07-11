from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system

class Lastz(AbstractMapper):
    def run(self):
        system("lastz %s %s --format=sam > %s" % (self.referenceFastaFile, self.readFastqFile, self.outputSamFile))
        
class LastzChain(Lastz):
    def run(self):
        Lastz.run(self)
        self.chainSamFile()
        
class LastzRealign(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile()