from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
from nanopore.analyses.utils import samToBamFile, samIterator
import os
import pysam

class QualiMap(AbstractAnalysis):
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        emptyQual = False
        for entry in samIterator(self.samFile):
            if entry.qual is None:
                emptyQual = True
        if emptyQual is False:
            localBamFile = os.path.join(self.getLocalTempDir(), "mapping.bam")
            localSortedBamFile = os.path.join(self.getLocalTempDir(), "mapping.sorted")
            samToBamFile(self.samFile, localBamFile)
            pysam.sort(localBamFile, localSortedBamFile)
            system("qualimap bamqc -bam %s -outdir %s" % (localSortedBamFile + ".bam", self.outputDir))
        self.finish()
