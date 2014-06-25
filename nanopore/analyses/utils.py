import pysam
from jobTree.src.bioio import reverseComplement, fastaRead

class AlignedPair:
    """Represents an aligned pair of positions.
    """
    def __init__(self, refPos, refSeq, readPos, isReversed, readSeq, pPair):
        assert refPos >= 0 and refPos < len(refSeq)
        self.refPos = refPos
        self.refSeq = refSeq
        assert readPos >= 0 and readPos < len(readSeq)
        self.readPos = readPos
        self.isReversed = isReversed
        self.readSeq = readSeq
        self.pPair = pPair #Pointer to the previous aligned pair
        
    def isMatch(self):
        return self.getRefBase().upper() == self.getReadBase().upper() and self.getRefBase().upper() in "ACTG"
    
    def isMismatch(self):
        return self.getRefBase().upper() != self.getReadBase().upper() and self.getRefBase().upper() in "ACTG" and self.getReadBase().upper() in "ACTG"
    
    def getRefBase(self):
        return self.refSeq[self.refPos]
    
    def getReadBase(self):
        if self.isReversed:
            return reverseComplement(self.readSeq[self.readPos]) 
        return self.readSeq[self.readPos]
    
    def getPrecedingReadInsertionLength(self):
        if self.pPair == None:
            return 0
        return self._indelLength(self.readPos, self.pPair.readPos)
    
    def getPrecedingReadDeletionLength(self):
        if self.pPair == None:
            return 0
        return self._indelLength(self.refPos, self.pPair.refPos)
    
    @staticmethod
    def _indelLength(pos, pPos):
        length = abs(pPos - pos) - 1
        assert length >= 0
        return length

    @staticmethod
    def iterator(alignedRead, refSeq, readSeq):
        """Generates aligned pairs from a pysam.AlignedRead object.
        """
        if alignedRead.cigar[0][0] == 5: #Translate the read position to the original coordinates by removing hard clipping
            readOffset = alignedRead.cigar[0][1]
        else:
            readOffset = 0
        if alignedRead.is_reverse: #SEQ is reverse complemented
            readOffset = -(len(readSeq) - 1 - readOffset)
        pPair = None
        for readPos, refPos in alignedRead.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedRead.pos and refPos < alignedRead.aend
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), alignedRead.is_reverse, readSeq, pPair)
                assert aP.getReadBase().upper() == alignedRead.seq[readPos].upper()
                pPair = aP
                yield aP

def samToBamFile(samInputFile, bamOutputFile):
    """Converts a sam file to a bam file (sorted)
    """
    #system("samtools view -Sb %s > %s" % (mapping.sam, mapping.bam))
    #system("samtools sort %s %s" % (mapping.bam, mapping.sorted))
    pass

