from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
from jobTree.src.bioio import system
from itertools import product

class UnmappedLengthDistributionAnalysis(AbstractUnmappedMetaAnalysis):
    """runs length distribution analysis on all mapped/unmapped per read Type
    as well as per reference"""
    def run(self):
        for readType in self.info.readTypes:
            unmapped = open(os.path.join(self.getLocalTempDir(), readType + "_unmapped"), "w")
            mapped = open(os.path.join(self.getLocalTempDir(), readType + "_mapped"), "w")
            for read in self.reads:
                if read.is_mapped:
                    mapped.write("{}\n".format(len(read.seq)))
                else:
                    unmapped.write("{}\n".format(len(read.seq)))
            system("Rscript nanopore/metaAnalyses/unmapped_mapped_distributions.R {} {} {} {}".format(os.path.join(self.getLocalTempDir(), readType + "_unmapped"), os.path.join(self.getLocalTempDir(), readType + "_mapped"), os.path.join(self.outputDir, readType + "_length_distribution.pdf"), readType))

        for reference in self.info.referenceFiles:
            unmapped = open(os.path.join(self.getLocalTempDir(), reference + "_unmapped"), "w")
            mapped = open(os.path.join(self.getLocalTempDir(), reference + "_mapped"), "w")
            for read in self.reads:
                if read.is_mapped:
                    mapped.write("{}\n".format(len(read.seq)))
                else:
                    unmapped.write("{}\n".format(len(read.seq)))
            system("Rscript nanopore/metaAnalyses/unmapped_mapped_distributions.R {} {} {} {}".format(os.path.join(self.getLocalTempDir(), reference + "_unmapped"), os.path.join(self.getLocalTempDir(), reference + "_mapped"), os.path.join(self.outputDir, reference + "_length_distribution.pdf"), reference))

        