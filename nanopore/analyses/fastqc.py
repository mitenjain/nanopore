from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
import os

class FastQC(AbstractAnalysis):
    def run(self):
        system("fastqc %s --outdir=%s" % (self.readFastqFile, self.outputDir))
