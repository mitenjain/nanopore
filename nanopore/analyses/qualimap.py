from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
from nanopore.analyses.utils import samToBamFile
import os
import pysam

class QualiMap(AbstractAnalysis):
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        localBamFile = os.path.join(self.getLocalTempDir(), "mapping.bam")
        localSortedBamFile = os.path.join(self.getLocalTempDir(), "mapping.sorted")
        samToBamFile(self.samFile, localBamFile)
        pysam.sort(localBamFile, localSortedBamFile)
        system("qualimap bamqc -bam %s -outdir %s" % (localSortedBamFile + ".bam", self.outputDir))
        self.finish()