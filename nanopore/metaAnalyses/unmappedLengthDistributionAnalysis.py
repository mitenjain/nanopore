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
            unmapped.close(); mapped.close()
            if os.path.getsize(os.path.join(self.getLocalTempDir(), readType + "_unmapped")) > 0 and os.path.getsize(os.path.join(self.getLocalTempDir(), readType + "_mapped")) > 0:
                system("Rscript nanopore/metaAnalyses/unmapped_mapped_distributions.R {} {} {} {}".format(os.path.join(self.getLocalTempDir(), readType + "_unmapped"), os.path.join(self.getLocalTempDir(), readType + "_mapped"), os.path.join(self.outputDir, readType + "_length_distribution.pdf"), readType))
            else:
                os.remove(os.path.join(self.getLocalTempDir(), readType + "_unmapped"))
                os.remove(os.path.join(self.getLocalTempDir(), readType + "_mapped"))
                
        for reference in self.info.referenceFiles:
            unmapped = open(os.path.join(self.getLocalTempDir(), os.path.basename(reference) + "_unmapped"), "w")
            mapped = open(os.path.join(self.getLocalTempDir(), os.path.basename(reference) + "_mapped"), "w")
            for read in self.reads:
                if read.is_mapped:
                    mapped.write("{}\n".format(len(read.seq)))
                else:
                    unmapped.write("{}\n".format(len(read.seq)))
            unmapped.close(); mapped.close()
            if os.path.getsize(os.path.join(self.getLocalTempDir(), readType + "_unmapped")) > 0 and os.path.getsize(os.path.join(self.getLocalTempDir(), readType + "_mapped")) > 0:
                system("Rscript nanopore/metaAnalyses/unmapped_mapped_distributions.R {} {} {} {}".format(os.path.join(self.getLocalTempDir(), os.path.basename(reference) + "_unmapped"), os.path.join(self.getLocalTempDir(), os.path.basename(reference) + "_mapped"), os.path.join(self.outputDir, os.path.basename(reference) + "_length_distribution.pdf"), os.path.basename(reference)))
            else:
                os.remove(os.path.join(self.getLocalTempDir(), readType + "_unmapped"))
                os.remove(os.path.join(self.getLocalTempDir(), readType + "_mapped"))