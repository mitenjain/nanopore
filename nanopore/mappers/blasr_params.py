from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class BlasrParams(AbstractMapper):
    def run(self):
        system("blasr %s %s -sdpTupleSize 8 -bestn 1 -nproc 8 -sam -out %s -m 0" % (self.readFastqFile, self.referenceFastaFile, self.outputSamFile))
