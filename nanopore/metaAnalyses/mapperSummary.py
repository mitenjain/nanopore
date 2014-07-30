from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os
import xml.etree.cElementTree as ET

class MapperSummary(AbstractMetaAnalysis):
    """Calculates meta-coverage across all the samples.
    """
    def run(self):
        fH = open(os.path.join(self.outputDir, "summary.tsv"), 'w')
        fH.write(",".join(["ReadFile", "ReferenceFile", "Mapper", "MedianReadCoverage","MedianReferenceCoverage","MedianIdentity","AveragePosteriorMatchProbability"]))
        for readFastqFile, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
            globalCoverageXML = ET.parse(os.path.join(resultsDir, "analysis_GlobalCoverage", "coverage_bestPerRead.xml")).getroot()
            alignmentUncertaintyXML = ET.parse(os.path.join(resultsDir, "analysis_AlignmentUncertainty", "alignmentUncertainty.xml")).getroot()
            indelsXML = ET.parse(os.path.join(resultsDir, "analysis_Indels", "indels.xml")).getroot()
            fH.write(",".join([readFastqFile, referenceFastaFile, mapper.__name__,
                               globalCoverageXML.attrib["medianreadCoverage"],
                               globalCoverageXML.attrib["medianreferenceCoverage"],
                               globalCoverageXML.attrib["medianidentity"],
                               alignmentUncertaintyXML.attrib["averagePosteriorMatchProbability"]]))
        fH.close()
        