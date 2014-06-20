from nanopore.src.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system

class Lastz(AbstractMapper):
    def run(self):
        system("lastz %s %s --format=sam > %s" % (self.referenceFastaFile, self.readFastaFile, self.outputSamFile))