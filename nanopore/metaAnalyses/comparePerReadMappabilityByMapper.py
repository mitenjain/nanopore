from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
from itertools import product
import re
from collections import OrderedDict as od

class ComparePerReadMappabilityByMapper(AbstractUnmappedMetaAnalysis):
    """Finds which base mappers mapped which reads"""
    def run(self):
        for readType in self.readTypes:
            sortedBaseMappers = sorted(self.baseMappers)
            outf = open(os.path.join(self.outputDir, readType + "_perReadMappability.tsv"), "w")
            outf.write("Read\tReadFastqFile\t"); outf.write("\t".join(sortedBaseMappers)); outf.write("\n")
            for read in self.reads:
                if read.readType == readType and read.is_mapped is True:
                    tmp = od([[x, 0] for x in sortedBaseMappers])
                    for mapper, reference in read.get_map_ref_pair():
                        baseMapper = re.findall("[A-Z][a-z]*", mapper)[0]
                        if tmp[baseMapper] == 0:
                            tmp[baseMapper] = 1
                    outf.write("\t".join([read.name, os.path.basename(read.readFastqFile)] + map(str, tmp.values()))); outf.write("\n")
            outf.close()
