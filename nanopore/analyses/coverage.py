from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import numpy
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, prettyXml, system

class ReadAlignmentCoverageCounter:
    """Counts coverage from a pairwise alignment.
    Global alignment means the entire reference and read sequences (trailing indels).
    """
    def __init__(self, readSeqName, readSeq, refSeqName, refSeq, alignedRead, globalAlignment=False):
        self.matches = 0
        self.mismatches = 0
        self.ns = 0
        self.totalReadInsertionLength = 0
        self.totalReadInsertions = 0
        self.totalReadDeletionLength = 0
        self.totalReadDeletions = 0
        self.readSeqName = readSeqName
        self.readSeq = readSeq
        self.refSeqName = refSeqName
        self.refSeq = refSeq
        self.globalAlignment = globalAlignment
        
        #Now process the read alignment
        totalReadInsertionLength, totalReadDeletionLength = 0, 0
        aP = None
        for aP in AlignedPair.iterator(alignedRead, self.refSeq, self.readSeq): 
            if aP.isMatch():
                self.matches += 1
            elif aP.isMismatch():
                self.mismatches += 1
            else:
                self.ns += 1
            if aP.getPrecedingReadInsertionLength(self.globalAlignment) > 0:
                self.totalReadInsertions += 1
                totalReadInsertionLength += aP.getPrecedingReadInsertionLength(self.globalAlignment)
            if aP.getPrecedingReadDeletionLength(self.globalAlignment) > 0:
                self.totalReadDeletions += 1
                totalReadDeletionLength += aP.getPrecedingReadDeletionLength(self.globalAlignment)
        if self.globalAlignment and aP != None: #If global alignment account for any trailing indels
            assert len(self.refSeq) - aP.refPos - 1 >= 0
            if len(self.refSeq) - aP.refPos - 1 > 0:
                self.totalReadDeletions += 1
                self.totalReadDeletionLength += len(self.refSeq) - aP.refPos - 1

            if alignedRead.is_reverse:
                aP.readPos >= 0
                if aP.readPos > 0:
                    self.totalReadInsertions += 1
                    totalReadInsertionLength += aP.readPos
            else:
                assert len(self.readSeq) - aP.readPos - 1 >= 0
                if len(self.readSeq) - aP.readPos - 1 > 0:
                    self.totalReadInsertions += 1
                    totalReadInsertionLength += len(self.readSeq) - aP.readPos - 1

        assert totalReadInsertionLength <= len(self.readSeq)
        assert totalReadDeletionLength <= len(self.refSeq)
        self.totalReadInsertionLength += totalReadInsertionLength
        self.totalReadDeletionLength += totalReadDeletionLength
            
    def readCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def referenceCoverage(self):
        return AbstractAnalysis.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def identity(self):
        return AbstractAnalysis.formatRatio(self.matches, self.matches + self.mismatches)
    
    def deletionsPerReadBase(self):
        return AbstractAnalysis.formatRatio(self.totalReadDeletions, len(self.readSeq))
    
    def insertionsPerReadBase(self):
        return AbstractAnalysis.formatRatio(self.totalReadInsertions, len(self.readSeq))
    
    def readLength(self):
        return len(self.readSeq)
    
    def getXML(self):
        return ET.Element("readAlignmentCoverage", { "refSeqName":self.refSeqName, 
                                                    "readSeqName":self.readSeqName, "readLength":str(self.readLength()),
                                                    "readCoverage":str(self.readCoverage()), 
                                                    "referenceCoverage":str(self.referenceCoverage()), 
                                                    "identity":str(self.identity()), 
                                                    "insertionsPerReadBase":str(self.insertionsPerReadBase()),
                                                    "deletionsPerReadBase":str(self.deletionsPerReadBase()) })

def getAggregateCoverageStats(readAlignmentCoverages, tagName, refSequences, readSequences, readsToReadAlignmentCoverages):
    """Calculates aggregate stats across a set of read alignments, plots distributions.
    """
    mappedReadLengths = [ len(readSequences[i]) for i in readSequences.keys() if i in readsToReadAlignmentCoverages ]
    unmappedReadLengths = [ len(readSequences[i]) for i in readSequences.keys() if i not in readsToReadAlignmentCoverages ]
    def stats(fnStringName):
        l = map(lambda x : getattr(x, fnStringName)(), readAlignmentCoverages)
        l.sort()
        return l[0], numpy.average(l), numpy.median(l), l[-1], " ".join(map(str, l))
    
    attribs = { "numberOfReadAlignments":str(len(readAlignmentCoverages)), 
                "numberOfReads":str(len(readSequences)), "numberOfReferenceSequences":str(len(refSequences)), 
                "numberOfMappedReads":str(len(mappedReadLengths)), "mappedReadLengths":" ".join(map(str, mappedReadLengths)),
                "numberOfUnmappedReads":str(len(unmappedReadLengths)), "unmappedReadLengths":" ".join(map(str, unmappedReadLengths)), }
    
    for fnStringName in "readCoverage", "referenceCoverage", "identity", "deletionsPerReadBase", "insertionsPerReadBase", "readLength":
        for attribName, value in zip([ "min" + fnStringName, "avg" + fnStringName, "median" + fnStringName, "max" + fnStringName, "distribution" + fnStringName ], list(stats(fnStringName))):
            attribs[attribName] = str(value)
            
    parentNode = ET.Element(tagName, attribs)
    for readAlignmentCoverage in readAlignmentCoverages:
        parentNode.append(readAlignmentCoverage.getXML())
    return parentNode

class LocalCoverage(AbstractAnalysis):
    """Calculates coverage, treating alignments as local alignments.
    """
    def run(self, globalAlignment=False):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        readsToReadCoverages = {}
        for aR in samIterator(sam): #Iterate on the sam lines
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = readSequences[aR.qname]
            readAlignmentCoverageCounter = ReadAlignmentCoverageCounter(aR.qname, readSeq, sam.getrname(aR.rname), refSeq, aR, globalAlignment)
            if aR.qname not in readsToReadCoverages:
                readsToReadCoverages[aR.qname] = []
            readsToReadCoverages[aR.qname].append(readAlignmentCoverageCounter)
        sam.close()
        #Write out the coverage info for differing subsets of the read alignments
        if len(readsToReadCoverages.values()) > 0:
            for readCoverages, outputName in [ (reduce(lambda x, y : x + y, readsToReadCoverages.values()), "coverage_all"), (map(lambda x : max(x, key=lambda y : y.readCoverage()), readsToReadCoverages.values()), "coverage_bestPerRead") ]:
                parentNode = getAggregateCoverageStats(readCoverages, outputName, refSequences, readSequences, readsToReadCoverages)
                open(os.path.join(self.outputDir, outputName + ".xml"), 'w').write(prettyXml(parentNode))
                #this is a ugly file format with each line being a different data type - column length is variable
                outf = open(os.path.join(self.outputDir, outputName + ".txt"), "w")
                outf.write("MappedReadLengths " + parentNode.get("mappedReadLengths") + "\n")
                outf.write("UnmappedReadLengths " + parentNode.get("unmappedReadLengths") + "\n")
                outf.write("ReadCoverage " + parentNode.get("distributionreadCoverage") + "\n")
                outf.write("ReadIdentity " + parentNode.get("distributionidentity") + "\n")
                outf.write("InsertionsPerBase " + parentNode.get("distributioninsertionsPerReadBase") + "\n")
                outf.write("DeletionsPerBase " + parentNode.get("distributiondeletionsPerReadBase") + "\n")
                outf.close()
                system("Rscript nanopore/analyses/coverage_plot.R {} {}".format(os.path.join(self.outputDir, outputName + ".txt"), os.path.join(self.outputDir, outputName + ".pdf")))
        self.finish()

class GlobalCoverage(LocalCoverage):
    def run(self):
        """Calculates coverage, treating alignments as global alignments.
        """
        LocalCoverage.run(self, globalAlignment=True)