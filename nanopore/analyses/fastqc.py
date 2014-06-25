from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
import os

class FastQC(AbstractMapper):
    def run(self):
        system("fastqc %s --outdir=%s" % (self.readFastQFile, self.outputDir))
