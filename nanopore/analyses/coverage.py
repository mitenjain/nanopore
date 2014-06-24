from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, prettyXml

class CoverageCounter:
    """Counts coverage from a pairwise alignment
    """
    def __init__(self):
        self.matches = 0
        self.mismatches = 0
        self.ns = 0
        self.totalReadInsertionLength = 0
        self.totalReadDeletionLength = 0
        
    def addReadAlignment(self, alignedRead, referenceSeq, readSeq):
        for aP in AlignedPair.iterator(alignedRead, refSeq, readSeq): 
            if aP.isMatch():
                self.matches += 1
            elif aP.isMismatch():
                self.mismatches += 1
            else:
                self.ns += 1
            self.totalReadInsertionLength += aP.getPrecedingReadInsertionLength()
            self.totalReadDeletionLength += aP.getPrecedingReadDeletionLength()
            
    def getReadCoverage(self):
        return self.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def getReferenceCoverage(self):
        return self.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def getAlignmentCoverage(self):
        return self.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    def getIdentity(self):
        return self.formatRatio(self.matches, self.matches + self.mismatches)
    
    def getReadIdentity(self):
        return self.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def getReferenceIdentity(self):
        return self.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def getAlignmentIdentity(self):
        return self.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    @staticmethod
    def _formatRatio(numerator, denominator):
        if denominator == 0:
            return float("nan")
        return float(numerator)/denominator
    
    def getXML(self):
        return ET.Element("coverage", { "readCoverage":str(self.getReadCoverage()), "referenceCoverage":str(self.getReferenceCoverage()), 
                                "alignmentCoverage":str(self.getAlignmentCoverage()), "identity":str(self.getIdentity()), 
                                "readIdentity":str(self.getReadIdentity()), "referenceIdentity":str(self.getReferenceIdentity()), 
                                "alignmentIdentity":str(self.getAlignmentIdentity()),
                                "matches":str(self.matches), "mismatches":str(self.mismatches), "ns":str(self.ns), 
                                "totalReadInsertionLength":str(self.totalReadInsertionLength),
                                "totalReadDeletionLength":str(self.totalReadDeletionLength) })

class Coverage(AbstractAnalysis):
    """Target calculates coverages.
    """
    def run(self):
        refSequences = dict(fastaRead(open(self.referenceFastaFile, 'r'))) #Hash of names to sequences
        readSequences = dict(fastaRead(open(self.readFastaFile, 'r'))) #Hash of names to sequences
        overallCoverageCounter = CoverageCounter() #Thing to store the overall coverage in
        readCoverages = []
        sam = pysam.Samfile(self.samFile, "r" )
        for aR in sam: #Iterate on the sam lines
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = readSequences[aR.qname]
            overallCoverageCounter.addReadAlignment(aR, refSeq, readSeq)
            readCoverages.append(CoverageCounter())
            readCoverages[-1].addReadAlignment(aR, refSeq, readSeq)   
        sam.close()
        #Write out the coverage info
        parentNode = overallCoverageCounter.getXML()
        for readCoverage in readCoverages:
            parentNode.add(readCoverage.getXML())
        open(os.path.join(self.outputDir, "coverages.xml"), 'w').write(prettyXml(parentNode))
        