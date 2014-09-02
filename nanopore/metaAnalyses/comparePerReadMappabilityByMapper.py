from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
from itertools import product
import re
from collections import Counter


class ComparePerReadMappabilityByMapper(AbstractUnmappedMetaAnalysis):
    """Finds which mappers mapped which reads; reports proportion of reads mapped by only one mapper"""
    def run(self):
        for readType in self.readTypes:
            mapper_ref_dict = {x: list() for x in self.baseMappers}
            for read in self.reads:
                if read.readType == readType and read.is_mapped and len(read.mapRefPairs) == 1:
                    for mapper, referenceFastaFile in read.get_map_ref_pair():
                        base_mapper = re.findall("[A-Z][a-z]*", mapper)[0]
                        mapper_ref_dict[base_mapper].add((read.mapRefPairs))

            outf = open(os.path.join(self.outputDir, readType + "_full_unique_results"), "w")
            outf2 = open(os.path.join(self.outputDir, readType + "_unique_read_counts"), "w")
            for base_mapper, singletons in mapper_ref_dict.iteritems():
                outf.write(base_mapper + " :"); outf.write("\t".join(singletons)); outf.write("\n")
                outf2.write(base_mapper + "\t" + str(len(singletons)) + "\n")
            outf.close(); outf2.close()

            #you want a map of read file and read to base mappers I.E. read1, file1: last=0, lastz=1, etc
            per_read_counts = list()
            for read in self.reads:
                if read.readType == readType and read.is_mapped:
                    for mapper, referenceFastaFile in read.get_map_ref_pair():
                        base_mapper_counter[(read.name, read.readFastqFile)] 