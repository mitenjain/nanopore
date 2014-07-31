from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os
import xml.etree.cElementTree as ET

class CoverageSummary(AbstractMetaAnalysis):
    """Calculates meta-coverage across all the samples.
    """
    def run(self):
        fH = open(os.path.join(self.outputDir, "summary.csv"), 'w')
        
        fH.write(",".join(["ReadFile", "ReferenceFile", "Mapper", "MedianReadCoverage","MedianReferenceCoverage","MedianIdentity","MedianDeletionsPerReadBase", "MedianInsertionsPerReadBase","AveragePosteriorMatchProbability"]) + "\n")
        
        for readFastqFile in self.readFastqFiles:
            for referenceFastaFile in self.referenceFastaFiles:
                ###Read coverage plot (x-axis: read coverage, y-axis: counts, series: mappers)
                    
                ###Unmapped reads plot (x-axis: mappers, y-axis: # of unmapped reads)
                    
                ###Median read coverage vs. avg read alignment uncertainty plot (x-axis: median read coverage, y-axis: avg. alignment uncertainty, series: mappers)
                    
                ###Median read coverage vs. median read identity plot (x-axis: median read coverage, y-axis: median identity, series: mappers)
                    
                ###Median read coverage vs. median number of indels per base (x-axis: median read coverage, y-axis: median number of insertions+deletions per base, series: mappers
                    
                for mapper in self.mappers:
                    analyses, resultsDir = self.experimentHash[(readFastqFile, referenceFastaFile, mapper)]
                    globalCoverageXML = ET.parse(os.path.join(resultsDir, "analysis_GlobalCoverage", "coverage_bestPerRead.xml")).getroot()
                    alignmentUncertaintyXML = ET.parse(os.path.join(resultsDir, "analysis_AlignmentUncertainty", "alignmentUncertainty.xml")).getroot()
                    fH.write(",".join([readFastqFile, referenceFastaFile, mapper.__name__,
                               globalCoverageXML.attrib["medianreadCoverage"], globalCoverageXML.attrib["medianreferenceCoverage"],
                               globalCoverageXML.attrib["medianidentity"], 
                               globalCoverageXML.attrib["mediandeletionsPerReadBase"],
                               globalCoverageXML.attrib["medianinsertionsPerReadBase"],
                               alignmentUncertaintyXML.attrib["averagePosteriorMatchProbability"]]) + "\n")
        fH.close()