from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class Blasr(AbstractMapper):
    def run(self):
        system("blasr %s %s -sam > %s" % (self.referenceFastaFile, self.readFastaFile, self.outputSamFile))