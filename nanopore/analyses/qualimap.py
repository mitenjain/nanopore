from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
from utils import samToBamFile
import os
import pysam

class QualiMap(AbstractAnalysis):
    def run(self):
        localBamFile = os.path.join(self.getLocalTempDir(), "mapping.bam")
        localSortedBamFile = os.path.join(self.getLocalTempDir(), "mapping.sorted")
        samToBamFile(self.samFile, localBamFile)
        pysam.sort(localBamFile, localSortedBamFile)
        system("qualimap bamqc -bam %s -outdir %s" % (localSortedBamFile + ".bam", self.outputDir))
