from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system
import re
from itertools import product, izip
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

    def by_mapper_readtype_reference(self):
        entry_map = {x : list() for x in product(self.baseMappers, self.readTypes, [os.path.basename(x) for x in self.referenceFastaFiles])}
        for entry in self.db:
            entry_map[(entry.base_mapper, entry.readType, entry.referenceFastaFile)].append(entry)
        for (base_mapper, readType, referenceFastaFile), entries in entry_map.iteritems():
            name = "_".join([base_mapper, readType, referenceFastaFile])
            self.write_file_analyze(entries, name)

    def by_mapper_readfile(self):
        entry_map = {x: list() for x in product(self.baseMappers, [os.path.basename(x[0]) for x in self.readFastqFiles])}
        for entry in self.db:
            entry_map[(entry.base_mapper, entry.readFastqFile)].append(entry)
        for (base_mapper, readFastqFile), entries in entry_map.iteritems():
            name = "_".join([base_mapper, readFastqFile])
            self.write_file_analyze(entries, name)

    def by_reference(self):
        entry_map = {os.path.basename(x): list() for x in self.referenceFastaFiles}
        for entry in self.db:
            entry_map[entry.referenceFastaFile].append(entry)
        for referenceFastaFile, entries in entry_map.iteritems():
            self.write_file_analyze(entries, referenceFastaFile, multiple_read_types=True)

    def write_file_analyze(self, entries, name, multiple_read_types=False):
        path = os.path.join(self.outputDir, name + ".csv")
        outf = open(path, "w")
        outf.write(",".join(["Name", "Mapper", "ReadType", "ReadFile", "ReferenceFile",  "AvgReadCoverage", "AvgReferenceCoverage", "AvgIdentity", "AvgMismatchesPerReadBase", "AvgDeletionsPerReadBase", "AvgInsertionsPerReadBase", "NumberOfMappedReads", "NumberOfUnmappedReads", "NumberOfReads"])); outf.write("\n")
        entries = sorted(entries, key = lambda x: (x.mapper, x.readType, x.readFastqFile))
        names = self.resolve_duplicate_rownames(entries, multiple_read_types)
        for entry, n in izip(entries, names):
            outf.write(",".join([n, entry.mapper, entry.readType, entry.readFastqFile, entry.referenceFastaFile,
                               entry.XML.attrib["avgreadCoverage"], entry.XML.attrib["avgreferenceCoverage"],
                               entry.XML.attrib["avgidentity"], entry.XML.attrib["avgmismatchesPerReadBase"], 
                               entry.XML.attrib["avgdeletionsPerReadBase"],
                               entry.XML.attrib["avginsertionsPerReadBase"],
                               entry.XML.attrib["numberOfMappedReads"],
                               entry.XML.attrib["numberOfUnmappedReads"],
                               entry.XML.attrib["numberOfReads"]]) + "\n")
        outf.close()
        path2 = os.path.join(self.outputDir, name + "_distribution.csv")
        outf = open(path2, "w")
        for entry, n in izip(entries, names):
            outf.write(",".join([n] + entry.XML.attrib["distributionidentity"].split())); outf.write("\n")
        outf.close()
        system("Rscript nanopore/metaAnalyses/coverageSummaryPlots.R {} {} {}".format(path, name, os.path.join(self.outputDir, name + "_summary_plots.pdf")))
        system("Rscript nanopore/metaAnalyses/coveragePlots.R {} {} {}".format(path2, name, os.path.join(self.outputDir, name + "_distribution.pdf")))


    def resolve_duplicate_rownames(self, entries, multiple_read_types=False):
        last_mapper = entries[0].mapper; count = 0; start = True
        names = list()
        for entry in entries:
            if multiple_read_types is True and entry.mapper + "_" + entry.readType == last_mapper:
                count += 1
                last_mapper = entry.mapper + "_" + entry.readType
                if start is not True:
                    names.append(entry.mapper + "_" + entry.readType + "." + str(count))
                else:
                    names.append(entry.mapper + "_" + entry.readType)
                    start = False
            elif multiple_read_types is True:
                last_mapper = entry.mapper + "_" + entry.readType
                names.append(entry.mapper + "_" + entry.readType)
                count = 1
            elif multiple_read_types is False and entry.mapper == last_mapper:
                count += 1
                last_mapper = entry.mapper
                if start is not True:
                    names.append(entry.mapper + "." + str(count))
                else:
                    names.append(entry.mapper)
                    start = False
            else:
                last_mapper = entry.mapper
                names.append(entry.mapper)
                count = 1
        return names


    def run(self):
        self.db = self.build_db()
        self.by_mapper_readtype_reference()
        self.by_mapper_readfile()
        self.by_reference()
