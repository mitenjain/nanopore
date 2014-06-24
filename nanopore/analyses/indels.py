from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, fastaRead, prettyXml

class Indels(AbstractAnalysis):
    """Calculates stats on indels.
    """
    def run(self):
        refSequences = dict(fastaRead(open(self.referenceFastaFile, 'r'))) #Hash of names to sequences
        readSequences = dict(fastaRead(open(self.readFastaFile, 'r'))) #Hash of names to sequences
        sM = SubstitutionMatrix() #The thing to store the counts in
        sam = pysam.Samfile(self.samFile, "r" )
        insertionsLengths = []
        deletionLengths = []
        blockLengths = []
        for aR in sam: #Iterate on the sam lines
            blockLength = 0
            for aP in AlignedPair.iterator(aR, refSequences[sam.getrname(aR.rname)], readSequences[aR.qname]):
                if aP.getPrecedingReadInsertionLength() > 0:
                    insertionLengths.append(aP.getPrecedingReadInsertionLength())
                if aP.getPrecedingReadDeletionLength() > 0:
                    deletionLengths.append(aP.getPrecedingReadDeletionLength())
                if aP.getPrecedingReadInsertionLength() > 0 or aP.getPrecedingReadDeletionLength() > 0:
                    assert blockLength > 0
                    blockLengths.append(blockLength)
                    blockLength = 1
                else:
                    blockLenth += 1
        sam.close()
        #Write out the substitution info
        open(os.path.join(self.outputDir, "indexl.xml"), 'w').write(prettyXml(sM.getXML()))
        