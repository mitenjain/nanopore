from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
import os

class QualiMap(AbstractAnalysis):
    def run(self):
	localBamFile = os.path.join(self.getLocalTempDir(), "mapping.bam")
	localSortedBamFile = os.path.join(self.getLocalTempDir(), "mapping.sorted")
	system("samtools view -Sb %s > %s" % (self.samFile, localBamFile))
	system("samtools sort %s %s" % (localBamFile, localSortedBamFile))
        system("qualimap bamqc -bam %s -outdir %s" % (localSortedBamFile + ".bam", self.outputDir))
