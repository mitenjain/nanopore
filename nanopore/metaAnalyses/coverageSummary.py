from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system

class CoverageSummary(AbstractMetaAnalysis):
    """Calculates meta-coverage across all the samples.
    """
    def run(self):
        fH = open(os.path.join(self.outputDir, "summary.csv"), 'w')
        
        fH.write(",".join(["ReadFile", "ReferenceFile", "Mapper", "AvgReadCoverage","AvgReferenceCoverage","AvgIdentity", "AvgMatchIdentity", "AvgDeletionsPerReadBase", "AvgInsertionsPerReadBase","AvgPosteriorMatchProbability", "UnmappedReadCount"]) + "\n")
        
        for readFastqFile in self.readFastqFiles:
            for referenceFastaFile in self.referenceFastaFiles:
                #tmp = open(os.path.join(self.outputDir, "tmp.csv"), "w")
                tmp = open(os.path.join(self.getLocalTempDir(), "tmp.csv"), "w")
                tmp.write(",".join(["Mapper", "AvgReadCoverage","AvgReferenceCoverage","AvgIdentity", "AvgMatchIdentity","AvgDeletionsPerReadBase", "AvgInsertionsPerReadBase","AvgPosteriorMatchProbability", "UnmappedReadCount", "NumberOfReads"]) + "\n")
                tmp_data = {}    
                for mapper in self.mappers:
                    analyses, resultsDir = self.experimentHash[(readFastqFile, referenceFastaFile, mapper)]
                    globalCoverageXML = ET.parse(os.path.join(resultsDir, "analysis_GlobalCoverage", "coverage_bestPerRead.xml")).getroot()
                    alignmentUncertaintyXML = ET.parse(os.path.join(resultsDir, "analysis_AlignmentUncertainty", "alignmentUncertainty.xml")).getroot()
                    if mapper.__name__ not in tmp_data:
                        tmp_data[mapper.__name__] = []
                    tmp_data[mapper.__name__].append(",".join(globalCoverageXML.attrib["distributionreadCoverage"].split()))
                    fH.write(",".join([readFastqFile, referenceFastaFile, mapper.__name__,
                               globalCoverageXML.attrib["avgreadCoverage"], globalCoverageXML.attrib["avgreferenceCoverage"],
                               globalCoverageXML.attrib["avgidentity"], globalCoverageXML.attrib["avgmatchIdentity"], 
                               globalCoverageXML.attrib["avgdeletionsPerReadBase"],
                               globalCoverageXML.attrib["avginsertionsPerReadBase"],
                               alignmentUncertaintyXML.attrib["averagePosteriorMatchProbability"],
                               globalCoverageXML.attrib["numberOfUnmappedReads"]]) + "\n")
                    #I realize this is ugly, sorry - better than trying to split it up in R however
                    tmp.write(",".join([mapper.__name__,
                               globalCoverageXML.attrib["avgreadCoverage"], globalCoverageXML.attrib["avgreferenceCoverage"],
                               globalCoverageXML.attrib["avgidentity"], globalCoverageXML.attrib["avgmatchIdentity"], 
                               globalCoverageXML.attrib["avgdeletionsPerReadBase"],
                               globalCoverageXML.attrib["avginsertionsPerReadBase"],
                               alignmentUncertaintyXML.attrib["averagePosteriorMatchProbability"],
                               globalCoverageXML.attrib["numberOfUnmappedReads"],
                               globalCoverageXML.attrib["numberOfReads"]]) + "\n")  
                tmp.close()
                filename = "_".join([readFastqFile.split("/")[-1],referenceFastaFile.split("/")[-1]])
                readname = readFastqFile.split("/")[-1]; refname = referenceFastaFile.split("/")[-1]
                #system("Rscript nanopore/metaAnalyses/coverageSummaryPlots.R {} {} {} {}".format(os.path.join(self.outputDir, "tmp.csv"), os.path.join(self.outputDir, filename), readname, refname))
                system("Rscript nanopore/metaAnalyses/coverageSummaryPlots.R {} {} {} {}".format(os.path.join(self.getLocalTempDir(), "tmp.csv"), os.path.join(self.outputDir, filename), readname, refname))
        fH.close()
        #tmp = open(os.path.join(self.outputDir, "tmp2.csv"), "w")
        tmp = open(os.path.join(self.getGlobalTempDir(), "tmp.csv"), "w")
        for mapper in tmp_data:
            tmp.write(",".join([mapper] + tmp_data[mapper])); tmp.write("\n")
            #Make version of the plot with just the mapper on it
            #tmp2 = open(os.path.join(self.getGlobalTempDir(), "tmp2.csv"), "w")
            #tmp2.write(",".join([mapper] + tmp_data[mapper])); tmp2.write("\n")
            #system("Rscript nanopore/metaAnalyses/coveragePlots.R {} {}".format(os.path.join(self.getGlobalTempDir(), "tmp2.csv"), os.path.join(self.outputDir, "coverage_summary_plots_%s.pdf" % mapper)))
            #tmp2.close()
        tmp.close()
        system("Rscript nanopore/metaAnalyses/coveragePlots.R {} {}".format(os.path.join(self.getGlobalTempDir(), "tmp.csv"), os.path.join(self.outputDir, "coverage_summary_plots.pdf")))
        #system("Rscript nanopore/metaAnalyses/coveragePlots.R {} {}".format(os.path.join(self.outputDir, "tmp2.csv"), os.path.join(self.outputDir, "coverage_summary_plots.pdf")))
