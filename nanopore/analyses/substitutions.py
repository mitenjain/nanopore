from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system
from itertools import product

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

    def getFreqs(self, refBase, bases):
        """
        Get list of relative frequencies for a refBase against all bases (passed as string)
        """
        freqs = list()
        for b in bases:
            freqs.append(self.getCount(refBase, b))
        if sum(freqs) == 0:
            return [ 0.0 ] * len(freqs)
        return [x / sum(freqs) for x in freqs]
    
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
    """Calculates stats on substitutions
    """
    def run(self, kmer=5):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sM = SubstitutionMatrix() #The thing to store the counts in
        sam = pysam.Samfile(self.samFile, "r" )
        for aR in samIterator(sam): #Iterate on the sam lines
            for aP in AlignedPair.iterator(aR, refSequences[sam.getrname(aR.rname)], readSequences[aR.qname]): #Walk through the matches mismatches:
                sM.addAlignedPair(aP.getRefBase(), aP.getReadBase())
        sam.close()
        #Write out the substitution info
        open(os.path.join(self.outputDir, "substitutions.xml"), 'w').write(prettyXml(sM.getXML()))
        bases = "ACGT"
        outf = open(os.path.join(self.outputDir, "subst.tsv"), "w")
        outf.write("A\tC\tG\tT\n")
        for x in bases:
            freqs = sM.getFreqs(x, bases)
            outf.write("{}\t{}\n".format(x, "\t".join(map(str,freqs)), "\n"))
        outf.close()
        analysis = self.outputDir.split("/")[-2].split("_")[-1] + "_Substitution_Levels"
        system("Rscript nanopore/analyses/substitution_plot.R {} {} {}".format(os.path.join(self.outputDir, "subst.tsv"), os.path.join(self.outputDir, "substitution_plot.pdf"), analysis))        
        self.finish()