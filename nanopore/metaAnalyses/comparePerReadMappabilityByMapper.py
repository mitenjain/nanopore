from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
from itertools import product
import re


class ComparePerReadMappabilityByMapper(AbstractUnmappedMetaAnalysis):
    """Finds which mappers mapped which reads; reports proportion of reads mapped by only one mapper"""
    def run(self):
        for readType in self.readTypes:
            mapper_ref_dict = {x: list() for x in self.baseMappers}
            for read in self.reads:
                if read.readType == readType and read.is_mapped and len(read.mapRefPairs) == 1:
                    for mapper, referenceFastaFile in read.get_map_ref_pair():
                        base_mapper = re.findall("[A-Z][a-z]*", mapper)[0]
                        mapper_ref_dict[base_mapper].add((read.name, read.readFastqFile))

        
            outf = open(os.path.join(self.outputDir, readType + "_full_unique_results"), "w")
            for base_mapper, singletons in mapper_ref_dict.iteritems():
                outf.write(base_mapper + " :"); outf.write("\t".join(singletons)); outf.write("\n")
            outf.close()

            outf = open(os.path.join(self.outputDir, readType + "_unique_read_counts"), "w")
            outf.write("\t".join(map(lambda x,y: x, str(len(y)), mapper_ref_dict.iteritems())) + "\n")
            outf.close()