from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, prettyXml, system

class CoverageCounter:
    """Counts coverage from a pairwise alignment
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
            self.totalReadDeletionLength +=  len(refSeq) - aP.refPos - 1
            if alignedRead.is_reverse:
                totalReadInsertionLength += aP.readPos
            else:
                assert len(readSeq) - aP.readPos - 1 >= 0
                totalReadInsertionLength += len(readSeq) - aP.readPos - 1
        assert totalReadInsertionLength <= len(readSeq)
        assert totalReadDeletionLength <= len(refSeq)
        self.totalReadInsertionLength += totalReadInsertionLength
        self.totalReadDeletionLength += totalReadDeletionLength
            
    def getReadCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def getReferenceCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def getAlignmentCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    def getIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches)
    
    def getReadIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def getReferenceIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def getAlignmentIdentity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength + self.totalReadDeletionLength)
    
    def getAlignedReferenceLength(self):
        return self.matches + self.mismatches + self.totalReadDeletionLength
    
    def getAlignedReadLength(self):
        return self.matches + self.mismatches + self.totalReadInsertionLength
    
    def getXML(self):
        return ET.Element("coverage", { "refSeqName":self.refSeqName, "readSeqName":self.readSeqName, "readCoverage":str(self.getReadCoverage()), "referenceCoverage":str(self.getReferenceCoverage()), 
                                "alignmentCoverage":str(self.getAlignmentCoverage()), "identity":str(self.getIdentity()), 
                                "readIdentity":str(self.getReadIdentity()), "referenceIdentity":str(self.getReferenceIdentity()), 
                                "alignmentIdentity":str(self.getAlignmentIdentity()),
                                "alignedReferenceLength":str(self.getAlignedReferenceLength()),
                                "alignedReadLength":str(self.getAlignedReadLength()),
                                "matches":str(self.matches), "mismatches":str(self.mismatches), "ns":str(self.ns), 
                                "totalReadInsertionLength":str(self.totalReadInsertionLength),
                                "totalReadDeletionLength":str(self.totalReadDeletionLength) })

class LocalCoverage(AbstractAnalysis):
    """Calculates coverage, treating alignments as local alignments.
    """
    def run(self, globalAlignment=False):
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        overallCoverageCounter = CoverageCounter("overall", "overall", globalAlignment=globalAlignment) #Thing to store the overall coverage in
        readCoverages = []
        sam = pysam.Samfile(self.samFile, "r" )
        for aR in samIterator(sam): #Iterate on the sam lines
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = readSequences[aR.qname]
            overallCoverageCounter.addReadAlignment(aR, refSeq, readSeq)
            readCoverages.append(CoverageCounter(aR.qname, sam.getrname(aR.rname), globalAlignment=globalAlignment))
            readCoverages[-1].addReadAlignment(aR, refSeq, readSeq)
        sam.close()
        #Write out the coverage info
        parentNode = overallCoverageCounter.getXML()
        for readCoverage in readCoverages:
            parentNode.append(readCoverage.getXML())
        open(os.path.join(self.outputDir, "coverages.xml"), 'w').write(prettyXml(parentNode))

        outf = open(os.path.join(self.outputDir, "coverages.tsv"), "w")
        outf.write("alignmentIdentity\talignmentCoverage\treadIdentity\treadCoverage\treferenceIdentity\treferenceCoverage\n")
        for readCoverage in readCoverages:
            outf.write("\t".join(map(str, [readCoverage.getAlignmentIdentity(), readCoverage.getAlignmentCoverage(), readCoverage.getReadIdentity(), readCoverage.getReadCoverage(), readCoverage.getReferenceIdentity(), readCoverage.getReferenceIdentity()])))
            outf.write("\n")
        outf.close()
        system("Rscript nanopore/analyses/coverage_plot.R {} {}".format(os.path.join(self.outputDir, "coverages.tsv"), os.path.join(self.outputDir, "coverage_histograms.pdf")))

class GlobalCoverage(LocalCoverage):
    def run(self):
        """Calculates coverage, treating alignments as global alignments.
        """
        LocalCoverage.run(self, globalAlignment=True)