from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system
import re
from collections import OrderedDict as od

class ComparePerReadMappabilityByMapper(AbstractUnmappedMetaAnalysis):
    """Finds which base mappers mapped which reads"""
    def run(self):
        for readType in self.readTypes:
            sortedBaseMappers = [x for x in sorted(self.baseMappers) if x != "Combined"]
            outf = open(os.path.join(self.outputDir, readType + "_perReadMappability.tsv"), "w")
            outf.write("Read\tReadFastqFile\t"); outf.write("\t".join(sortedBaseMappers)); outf.write("\n")
            for read in self.reads:
                if read.readType == readType:
                    tmp = od([[x, 0] for x in sortedBaseMappers])
                    if read.is_mapped is True:
                        for mapper, reference in read.get_map_ref_pair():
                            baseMapper = re.findall("[A-Z][a-z]*", mapper)[0]
                            #hacky way to avoid including 'combined' analysis
                            if baseMapper != "Combined" and tmp[baseMapper] == 0:
                                tmp[baseMapper] = 1
                    outf.write("\t".join([read.name, os.path.basename(read.readFastqFile)] + map(str, tmp.values()))); outf.write("\n")
            outf.close()
            system("Rscript nanopore/metaAnalyses/vennDiagram.R {} {}".format(os.path.join(self.outputDir, readType + "_perReadMappability.tsv"), os.path.join(self.outputDir, readType + "_perReadMappabilityVennDiagram.pdf")))