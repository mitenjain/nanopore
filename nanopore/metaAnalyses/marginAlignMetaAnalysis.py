from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead, fastaRead
from nanopore.analyses.utils import samIterator

class MarginAlignMetaAnalysis(AbstractMetaAnalysis):
    def run(self):
        for readFastqFile, readType in self.readFastqFiles:
            for referenceFastaFile in self.referenceFastaFiles:
                fH = open(os.path.join(self.outputDir, "marginAlign_%s_%s_%s.txt" % (os.path.split(readFastqFile)[-1], readType, os.path.split(referenceFastaFile)[-1])), 'w')
                fH.write("\t".join(["mapper", "caller", "totalHeldOut", "totalNonHeldOut", "totalHeldOutCallsTrue", "totalHeldOutCallsFalseAndReference", "totalHeldOutCallsFalseAndNonReference", "totalHeldOutNotCalled", "totalNonHeldOutCallsTrue", "totalNonHeldOutCallsFalse", "totalNonHeldOutNotCalled" ]) + "\n")
                for mapper in self.mappers:
                    analyses, resultsDir = self.experimentHash[((readFastqFile, readType), referenceFastaFile, mapper)]
                    node = ET.parse(os.path.join(resultsDir, "analysis_MarginAlignSnpCaller", "marginaliseConsensus.xml")).getroot()
                    for c in node:
                        fH.write("\t".join([mapper.__name__, c.tag, c.attrib["totalHeldOut"], c.attrib["totalNonHeldOut"], 
                                            c.attrib["totalHeldOutCallsTrue"], c.attrib["totalHeldOutCallsFalseAndReference"], 
                                            c.attrib["totalHeldOutCallsFalseAndNonReference"], c.attrib["totalHeldOutNotCalled"], 
                                            c.attrib["totalNonHeldOutCallsTrue"], c.attrib["totalNonHeldOutCallsFalse"], c.attrib["totalNonHeldOutNotCalled"] ]) + "\n")
                