from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system
from itertools import product
import pysam


class ReadTypeIntersectionsByMapper(AbstractMetaAnalysis):
    """Calculates the number of reads that mapped as twoD that did not map as template+complement by mapper/parameter set 
    """
    def run(self):
        combinedResults = dict()
        orderedReadTypes = sorted(self.readTypes)
        for mapper, (readFastqFile, readType), referenceFastaFile in product(self.mappers, self.readFastqFiles, self.referenceFastaFiles):
            readSets = dict()
            path = os.path.join("./output", "analysis_" + readType, "experiment_" + readFastqFile + "." + referenceFastaFile + "_" + mapper.__name__, "mapping.sam")
            if os.path.exists(path) is False:
                path = os.path.join("./tests/output", "analysis_" + readType, "experiment_" + readFastqFile + "." + referenceFastaFile + "_" + mapper.__name__, "mapping.sam")
            if os.path.exists(path) is False:
                break
            readSets[readType] = {x.qname for x in pysam.Samfile(path) if not x.is_unmapped}
            combinedResults[(mapper, readFastqFile, referenceFastaFile)] = readSets
        if sum([len(x) for x in combinedResults.values()]) > 0:
            outf = open(os.path.join(self.outputDir, "intersection_metrics.tsv"), "w")
            outf.write("Mapper\tReadFastqFile\tReferenceFastaFile\t{0}_Intersection_{1}\t{0}_Intersection_{2}\t{1}_Intersection_{2}\n".format(*orderedReadTypes))
            for (mapper, readFastqFile, referenceFastaFile), readSets in sorted(combinedResults.iteritems(), key = lambda x: x[0][0]):
                is1 = readSets[orderedReadTypes[0]].intersection(readSets[orderedReadTypes[1]])
                is2 = readSets[orderedReadTypes[0]].intersection(readSets[orderedReadTypes[2]])
                is3 = readSets[orderedReadTypes[1]].intersection(readSets[orderdReadTypes[2]])
                outf.write("\t".join(mapper.__name__, readFastqFile, referenceFastaFile, is1, is2, is3))
            outf.close()