from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system
import re
from itertools import product
from collections import Counter

class Entry(object):
    def __init__(self, readType, readFastqFile, referenceFastaFile, mapper, XML):
        self.readType = readType
        self.readFastqFile = os.path.basename(readFastqFile)
        self.referenceFastaFile = os.path.basename(referenceFastaFile)
        self.mapper = mapper.__name__
        self.base_mapper = re.findall("[A-Z][a-z]*", mapper.__name__)[0]
        self.XML = XML

class CoverageSummary(AbstractMetaAnalysis):
    """Calculates meta-coverage across all the samples.
    Includes analysis per mapper, per readType+mapper, per reference, and combined
    """
    def build_db(self):
        """
        Builds db of coverage results
        """
        db = list()
        for readFastqFile, readType in self.readFastqFiles:
            for referenceFastaFile in self.referenceFastaFiles:
                for mapper in self.mappers:
                    analyses, resultsDir = self.experimentHash[(readFastqFile, readType), referenceFastaFile, mapper]
                    if os.path.exists(os.path.join(resultsDir, "analysis_GlobalCoverage", "coverage_bestPerRead.xml")):
                        globalCoverageXML = ET.parse(os.path.join(resultsDir, "analysis_GlobalCoverage", "coverage_bestPerRead.xml")).getroot()
                        db.append(Entry(readType, readFastqFile, referenceFastaFile, mapper, globalCoverageXML))
        return db      

    def by_mapper_readtype(self):
        base_mappers = {re.findall("[A-Z][a-z]*", mapper.__name__)[0] for mapper in self.mappers}
        entry_map = {x : list() for x in product(base_mappers, self.readTypes)}
        for entry in self.db:
            entry_map[(entry.base_mapper, entry.readType)].append(entry)
        for (base_mapper, readType), entries in entry_map.iteritems():
            name = "_".join([os.path.basename(base_mapper), readType])
            self.write_file_analyze(entries, name)
        return entry_map

    #def by_reference(self):
    #    entry_map = {x: list() for x in self.referenceFastaFiles}
    #    for entry in self.db:
    #        entry_map[entry.referenceFastaFile] = entry
    #    return entry_map

    def write_file_analyze(self, entries, name):
        path = os.path.join(self.outputDir, name + ".tsv")
        outf = open(path, "w")
        outf.write(",".join(["Mapper", "ReadFile", "ReferenceFile",  "AvgReadCoverage", "AvgReferenceCoverage", "AvgIdentity", "AvgMatchIdentity", "AvgDeletionsPerReadBase", "AvgInsertionsPerReadBase", "NumberOfMappedReads", "NumberOfUnmappedReads", "NumberOfReads"])); outf.write("\n")
        entries = self.resolve_duplicate_rownames(entries)
        for entry in entries:
            outf.write(",".join([entry.mapper, entry.readFastqFile, entry.referenceFastaFile,
                               entry.XML.attrib["avgreadCoverage"], entry.XML.attrib["avgreferenceCoverage"],
                               entry.XML.attrib["avgidentity"], entry.XML.attrib["avgmatchIdentity"], 
                               entry.XML.attrib["avgdeletionsPerReadBase"],
                               entry.XML.attrib["avginsertionsPerReadBase"],
                               entry.XML.attrib["numberOfMappedReads"],
                               entry.XML.attrib["numberOfUnmappedReads"],
                               entry.XML.attrib["numberOfReads"]]) + "\n")
        outf.close()
        path2 = os.path.join(self.outputDir, name + "_distribution.tsv")
        outf = open(path2, "w")
        for entry in entries:
            outf.write(",".join([entry.mapper] + entry.XML.attrib["distributionidentity"].split())); outf.write("\n")
        outf.close()
        system("Rscript nanopore/metaAnalyses/coverageSummaryPlots.R {} {} {}".format(path, name, os.path.join(self.outputDir, name)))
        system("Rscript nanopore/metaAnalyses/coveragePlots.R {} {} {}".format(path2, name, os.path.join(self.outputDir, name + "_distribution")))


    def resolve_duplicate_rownames(self, entries):
        mappers = Counter()
        for entry in entries:
            if entry.mapper not in mappers:
                mappers[entry.mapper] += 1
            else:
                entry.mapper = entry.mapper + "." + str(mappers[entry.mapper])
                mappers[entry.mapper] += 1
        return entries


    def run(self):
        self.db = self.build_db()
        self.by_mapper_readtype()

