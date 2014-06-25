from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
import os

class QualiMap(AbstractAnalysis):
    def run(self):
    	mapping.sam = self.samFile
	system("samtools view -Sb %s > %s" % (mapping.sam, mapping.bam))
	system("samtools sort %s %s" % (mapping.bam, mapping.sorted))
        system("qualimap/qualimap bamqc -bam -outdir %s %s" % (mapping.sorted.bam, self.outputDir))

