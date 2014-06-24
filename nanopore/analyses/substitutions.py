from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, prettyXml

class SubstitutionMatrix():
    """Represents a nucleotide substitution matrix. Also allows 
    for recording matches against Ns.
    """
    def __init__(self):
        self.matrix = [0.0]*25 #Includes alignments between wildcard characters.
    
    def addAlignedPair(self, refBase, readBase):
        self.matrix[self._index(refBase) * 5 + self._index(readBase)] += 1
    
    def getCount(self, refBase, readBase):
        return self.matrix[self._index(refBase) * 5 + self._index(readBase)]
    
    def getXML(self):
        def _identity(matches, mismatches):
            if matches + mismatches == 0:
                return "NaN"
            return matches/(mismatches+matches)
        matches = sum([ self.getCount(base, base) for base in "ACTG" ])
        mismatches = sum([ sum([ self.getCount(refBase, readBase) for readBase in "ACTG" if readBase != refBase ]) for refBase in "ACTG" ])
        node = ET.Element("substitutions", { "matches":str(matches), "mismatches":str(mismatches), "identity":str(_identity(matches, mismatches)) })
        overallMatches = 0
        overallMismatches = 0
        for refBase in "ACGTN":
            matches = self.getCount(refBase, refBase)
            mismatches = sum([ self.getCount(refBase, readBase) for readBase in "ACTG" if readBase != refBase ])
            baseNode = ET.SubElement(node, refBase, { "matches":str(matches), "mismatches":str(mismatches), "identity":str(_identity(matches, mismatches)) })
            for readBase in "ACGTN":
                ET.SubElement(baseNode, readBase, { "count":str(self.getCount(refBase, readBase)) })
        return node
    
    @staticmethod
    def _index(base):
        base = base.upper()
        if base not in "ACGT":
            return 4
        return { 'A':0, 'C':1, 'G':2, 'T':3 }[base]

class Substitutions(AbstractAnalysis):
    """This is just my first test analysis target
    """
    def run(self):
        refSequences = dict(fastaRead(open(self.referenceFastaFile, 'r'))) #Hash of names to sequences
        readSequences = dict(fastaRead(open(self.readFastaFile, 'r'))) #Hash of names to sequences
        dM = DistanceMatrix() #The thing to store the counts in
        sam = pysam.Samfile(self.samFile, "r" )
        for aR in sam: #Iterate on the sam lines
            for aP in alignedPairIt(aR, refSequences[sam.getrname(aR.rname)], readSequences[aR.qname]): #Walk through the matches mismatches:
                dM.recordAlignedPair(aP.getRefBase(), aP.getReadBase())
        sam.close()
        #Write out the substitution info
        open(os.path.join(self.outputDir, "substitutions.xml"), 'w').write(prettyXml(dM.getXML()))
        