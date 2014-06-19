from nanopore.src.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system

class Lastz(AbstractMapper):
    def run(self):
        system("cactus_lastz %s %s --format=cigar > %s" % (self.referenceFastaFile, self.readFastaFile, self.outputSamFile))