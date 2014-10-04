from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead, fastaRead
from nanopore.analyses.utils import samIterator
from itertools import product
import numpy

class MarginAlignMetaAnalysis(AbstractMetaAnalysis):
    def run(self):
        readTypes = set([ readType for readFastqFile, readType in self.readFastqFiles ])
        coverageLevels = set()
        hash = {}
        variantCallingAlgorithms = set()
        proportionsHeldOut = set()
        for referenceFastaFile in self.referenceFastaFiles:
            for readType in readTypes:
                for readFastqFile, readFileReadType in self.readFastqFiles:
                    if readFileReadType == readType:
                        for mapper in self.mappers:
                            analyses, resultsDir = self.experimentHash[((readFastqFile, readType), referenceFastaFile, mapper)]
                            try:
                                node = ET.parse(os.path.join(resultsDir, "analysis_MarginAlignSnpCaller", "marginaliseConsensus.xml")).getroot()
                            except:
                                self.logToMaster("Parsing the file following failed: %s" % os.path.join(resultsDir, "analysis_MarginAlignSnpCaller", "marginaliseConsensus.xml"))
                                continue
                            for c in node:
                                coverage = int(c.attrib["coverage"])
                                coverageLevels.add(coverage)
                                proportionHeldOut = float(c.attrib["totalHeldOut"]) / (float(c.attrib["totalHeldOut"]) + float(c.attrib["totalNonHeldOut"]))
                                key = (readType, mapper.__name__, c.tag, proportionHeldOut, referenceFastaFile)
                                variantCallingAlgorithms.add(c.tag)
                                proportionsHeldOut.add(proportionHeldOut)
                                if key not in hash:
                                    hash[key] = {}
                                if coverage not in hash[key]:
                                    hash[key][coverage] = []
                                hash[key][coverage].append(c)
                                
        fH = open(os.path.join(self.outputDir, "marginAlignAll.txt"), 'w')
        fH.write("\t".join(["readType", "mapper", "caller", 
                            "%heldOut", "coverage", 
                            "fScoreMin", "fScoreMedian", "fScoreMax",
                            "recallMin", "recallMedian", "recallMax",
                            "precisionMin", "precisionMedian", "precisionMax", 
                            "%notCalledMin", "%notCalledMedian", "%notCalledMax",
                            "actualCoverageMin", "actualCoverageMedian", "actualCoverageMax"]) + "\n")
        
        fH2 = open(os.path.join(self.outputDir, "marginAlignSquares.txt"), 'w')
        coverageLevels = list(coverageLevels)
        coverageLevels.sort()
        fH2.write("\t".join(["readType", "mapper", "caller", 
                            "%heldOut",
                           "\t".join([ ("min_recall_coverage_%s\t" % coverage) + ("avg_recall_coverage_%s\t" % coverage) + ("max_recall_coverage_%s" % coverage) for coverage in coverageLevels]),
                           "\t".join([ ("min_precision_coverage_%s\t" % coverage) + ("avg_precision_coverage_%s\t" % coverage) + ("max_precision_coverage_%s" % coverage) for coverage in coverageLevels]),
                           "\t".join([ ("min_fscore_coverage_%s\t" % coverage) + ("avg_fscore_coverage_%s\t" % coverage) + ("max_fscore_coverage_%s" % coverage) for coverage in coverageLevels]), "\n" ]))
        
        keys = hash.keys()
        keys.sort()
        
        rocCurvesHash = {}
        
        for readType, mapper, algorithm, proportionHeldOut, referenceFastaFile in keys:
            nodes = hash[(readType, mapper, algorithm, proportionHeldOut, referenceFastaFile)]
            
            recall = lambda c : float(c.attrib["recall"]) #float(c.attrib["totalTruePositives"])/float(c.attrib["totalHeldOut"]) if float(c.attrib["totalHeldOut"]) != 0 else 0
            precision = lambda c : float(c.attrib["precision"]) # float(c.attrib["totalTruePositives"])/(float(c.attrib["totalTruePositives"]) + float(c.attrib["totalFalsePositives"])) if float(c.attrib["totalTruePositives"]) + float(c.attrib["totalFalsePositives"]) != 0 else 0
            fScore = lambda c : 2 * precision(c) * recall(c) / (precision(c) + recall(c)) if precision(c) + recall(c) > 0 else 0
            notCalled = lambda c : float(c.attrib["totalNoCalls"]) / (float(c.attrib["totalHeldOut"]) + float(c.attrib["totalNonHeldOut"]))
            actualCoverage = lambda c : float(c.attrib["actualCoverage"])
            
            for coverage in coverageLevels:
                def r(f):
                    i = map(f, nodes[coverage])
                    return "\t".join(map(str, (min(i), numpy.median(i), max(i))))
                fH.write("\t".join([readType, mapper, algorithm, str(proportionHeldOut), str(coverage), 
                                   r(fScore), r(recall), r(precision), r(notCalled), r(actualCoverage)]) + "\n")
            
            fH2.write("\t".join([readType, mapper, algorithm, str(proportionHeldOut)]) + "\t")
            f2 = lambda f, end : fH2.write("\t".join(map(lambda coverage : str(min(map(f, nodes[coverage]))) + "\t" + str(numpy.average(map(f, nodes[coverage]))) + "\t" + str(max(map(f, nodes[coverage]))), coverageLevels)) + end)
            f2(recall, "\t")
            f2(precision, "\t")
            f2(fScore, "\n")
            
            #Make ROC curves
            for coverage in coverageLevels:
                #Get the median true positive / median false positives
                recallByProbability = map(lambda c : map(float, c.attrib["recallByProbability"].split()), nodes[coverage])
                precisionByProbability = map(lambda c : map(float, c.attrib["precisionByProbability"].split()), nodes[coverage])
                def merge(curves, fn):
                    return map(lambda i : fn(map(lambda curve : curve[i], curves)), range(len(curves[0])))
                avgRecallByProbability = merge(recallByProbability, numpy.average)
                avgPrecisionByProbability = merge(precisionByProbability, numpy.average)
                rocCurvesHash[(readType, mapper, algorithm, proportionHeldOut, coverage)] = (avgPrecisionByProbability, avgRecallByProbability)
        
        
        ####Ian todo ###
        
        #Place to create ROC / precision/recall plots
        variantCallingAlgorithms = list(variantCallingAlgorithms)
        variantCallingAlgorithms.sort()
        proportionsHeldOut = list(proportionsHeldOut)
        proportionsHeldOut.sort()
        for readType, mapper in product(readTypes, self.mappers):
            outf = open(os.path.join(self.outputDir, readType + "_" + mapper.__name__ + ".tsv"), "w")
            #Make grid plot for each combination of readType/mapper
            #Grid dimensions would be variant calling algorithms x proportion held out
            #On each plot we should show the roc curve (use falsePositiveRatesByProbability vs. truePositiveRatesByProbability) for the different coverages.
            for algorithm in variantCallingAlgorithms:
                for proportionHeldOut in proportionsHeldOut:
                    for coverage in coverageLevels:
                        avgPrecisionByProbability, avgRecallByProbability = rocCurvesHash[(readType, mapper.__name__, algorithm, proportionHeldOut, coverage)]
                        outf.write("FPR\t{0}\t{1}\t{2}\t{3}\nTPR\t{0}\t{1}\t{2}\t{4}\n".format(str(algorithm), str(proportionHeldOut), str(coverage), "\t".join(map(str,recallByProbability)), "\t".join(map(str,precisionByProbability))))
            outf.close()
            if not os.path.exists(os.path.join(self.outputDir, readType + "_" + mapper.__name__)):
                os.mkdir(os.path.join(self.outputDir, readType + "_" + mapper.__name__))
            system("Rscript nanopore/metaAnalyses/ROC_marginAlign.R {} {} {}".format(os.path.join(self.outputDir, readType + "_" + mapper.__name__ + ".tsv"), os.path.join(self.outputDir, readType + "_" + mapper.__name__) + "/", "_ROC_curves.pdf"))
