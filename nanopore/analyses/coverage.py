from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import numpy
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, prettyXml, system

class CoverageCounter:
    """Counts coverage from a pairwise alignment.
    Global alignment means the entire reference and read sequences (trailing indels).
    """
    def __init__(self, readSeqName, refSeqName, globalAlignment=False):
        self.matches = 0
        self.mismatches = 0
        self.ns = 0
        self.totalReadInsertionLength = 0
        self.totalReadDeletionLength = 0
        self.readSeqName = readSeqName
        self.refSeqName = refSeqName
        self.globalAlignment = globalAlignment
        
    def addReadAlignment(self, alignedRead, refSeq, readSeq):
        totalReadInsertionLength, totalReadDeletionLength = 0, 0
        for aP in AlignedPair.iterator(alignedRead, refSeq, readSeq): 
            if aP.isMatch():
                self.matches += 1
            elif aP.isMismatch():
                self.mismatches += 1
            else:
                self.ns += 1
            totalReadInsertionLength += aP.getPrecedingReadInsertionLength(self.globalAlignment)
            totalReadDeletionLength += aP.getPrecedingReadDeletionLength(self.globalAlignment)
        if self.globalAlignment: #If global alignment account for any trailing indels
            assert len(refSeq) - aP.refPos - 1 >= 0
            self.totalReadDeletionLength += len(refSeq) - aP.refPos - 1
            if alignedRead.is_reverse:
                totalReadInsertionLength += aP.readPos
            else:
                assert len(readSeq) - aP.readPos - 1 >= 0
                totalReadInsertionLength += len(readSeq) - aP.readPos - 1
        assert totalReadInsertionLength <= len(readSeq)
        assert totalReadDeletionLength <= len(refSeq)
        self.totalReadInsertionLength += totalReadInsertionLength
        self.totalReadDeletionLength += totalReadDeletionLength
            
    def readCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def referenceCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def alignmentCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    def identity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches)
    
    def readIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def referenceIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def alignmentIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    def alignedReferenceLength(self):
        return self.matches + self.mismatches + self.totalReadDeletionLength
    
    def alignedReadLength(self):
        return self.matches + self.mismatches + self.totalReadInsertionLength
    
    def getXML(self):
        return ET.Element("coverage", { "refSeqName":self.refSeqName, "readSeqName":self.readSeqName, "readCoverage":str(self.readCoverage()), "referenceCoverage":str(self.referenceCoverage()), 
                                "alignmentCoverage":str(self.alignmentCoverage()), "identity":str(self.identity()), 
                                "readIdentity":str(self.readIdentity()), "referenceIdentity":str(self.referenceIdentity()), 
                                "alignmentIdentity":str(self.alignmentIdentity()),
                                "alignedReferenceLength":str(self.alignedReferenceLength()),
                                "alignedReadLength":str(self.alignedReadLength()),
                                "matches":str(self.matches), "mismatches":str(self.mismatches), "ns":str(self.ns), 
                                "totalReadInsertionLength":str(self.totalReadInsertionLength),
                                "totalReadDeletionLength":str(self.totalReadDeletionLength) })

def getAggregateCoverageStats(readCoverages, tagName, refSequences, readSequences, readsToReadCoverages):
    """Calculates aggregate stats across a set of read alignments
    """
    def stats(fnStringName):
        l = map(lambda x : getattr(x, fnStringName)(), readCoverages)
        l.sort()
        return l[0], numpy.average(l), numpy.median(l), l[-1], " ".join(map(str, l))
    attribs = { "numberOfReadAlignments":str(len(readCoverages)), "numberOfReads":str(len(readSequences)), 
               "numberOfReferenceSequences":str(len(refSequences)), "numberOfReadsWithoutAnyAlignment":len(readsToReadCoverages) }
    for fnStringName in "readCoverage", "referenceCoverage", "alignmentCoverage", "identity", "readIdentity", "referenceIdentity", "alignmentIdentity":
        for attribName, value in zip([ "min" + fnStringName, "avg" + fnStringName, "median" + fnStringName, "max" + fnStringName, "distribution" + fnStringName ], list(stats(fnStringName))):
            attribs[attribName] = str(value)
    parentNode = ET.Element(tagName, attribs)
    for readCoverage in readCoverages:
            parentNode.append(readCoverage.getXML())
    return parentNode

class LocalCoverage(AbstractAnalysis):
    """Calculates coverage, treating alignments as local alignments.
    """
    def run(self, globalAlignment=False):
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        readsToReadCoverages = {}
        for aR in samIterator(sam): #Iterate on the sam lines
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = readSequences[aR.qname]
            coverageCounter = CoverageCounter(aR.qname, sam.getrname(aR.rname), globalAlignment=globalAlignment)
            coverageCounter.addReadAlignment(aR, refSeq, readSeq)
            if aR.qname not in readsToReadCoverages:
                readsToReadCoverages[aR.qname] = []
            readsToReadCoverages[aR.qname].append(coverageCounter)
        sam.close()
        #Write out the coverage info for differing subsets of the read alignments
        for readCoverages, outputName in [ (reduce(lambda x, y : x + y, readsToReadCoverages.values()), "coverage_all"), (map(lambda x : max(x, key=lambda y : y.readCoverage()), readsToReadCoverages.values()), "coverage_bestPerRead") ]:
            parentNode = getAggregateCoverageStats(readCoverages, outputName, refSequences, readSequences, readsToReadCoverages)
            open(os.path.join(self.outputDir, outputName + ".xml"), 'w').write(prettyXml(parentNode))
    
            if len(readCoverages) > 0:
                outf = open(os.path.join(self.outputDir, outputName + ".tsv"), "w")
                outf.write("alignmentIdentity\talignmentCoverage\treadIdentity\treadCoverage\treferenceIdentity\treferenceCoverage\n")
                for readCoverage in readCoverages:
                    outf.write("\t".join(map(str, [readCoverage.alignmentIdentity(), readCoverage.alignmentCoverage(), readCoverage.readIdentity(), readCoverage.readCoverage(), readCoverage.referenceIdentity(), readCoverage.referenceIdentity()])))
                    outf.write("\n")
                outf.close()
                system("Rscript nanopore/analyses/coverage_plot.R {} {}".format(os.path.join(self.outputDir, outputName + ".tsv"), os.path.join(self.outputDir, outputName + ".pdf")))

class GlobalCoverage(LocalCoverage):
    def run(self):
        """Calculates coverage, treating alignments as global alignments.
        """
        LocalCoverage.run(self, globalAlignment=True)